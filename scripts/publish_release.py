import argparse
from dataclasses import dataclass
import hashlib
from http.client import HTTPException, HTTPSConnection
import json
import os
from pathlib import Path
import re
from urllib.parse import quote, urlsplit


PRODUCT = "rsomics-fastq-preprocess"
VERSION = "0.2.0"
TAG = f"v{VERSION}"
REPOSITORY = f"omics-rust/{PRODUCT}"
TAG_REF = f"refs/tags/{TAG}"
API_ROOT = f"/repos/{REPOSITORY}"
UPLOAD_HOST = "uploads.github.com"
FILES = {f"{PRODUCT}-{VERSION}-{target}.tar.gz" for target in (
    "linux-x86_64", "linux-aarch64", "macos-x86_64", "macos-aarch64",
)} | {f"{PRODUCT}-{VERSION}.crate", "SHA256SUMS"}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(body):
    return hashlib.sha256(body).hexdigest()


@dataclass(frozen=True)
class TagRun:
    head: str
    run_id: int
    attempt: int

    @classmethod
    def read(cls, env, event):
        require(env.get("GITHUB_EVENT_NAME") == "push" and env.get("GITHUB_REF") == TAG_REF
                and env.get("GITHUB_REPOSITORY") == REPOSITORY and event.get("ref") == TAG_REF
                and event.get("created") is True and event.get("deleted") is False,
                f"only a new {TAG} tag push can publish this release")
        head = env["GITHUB_SHA"]
        run_id, attempt = int(env["GITHUB_RUN_ID"]), int(env["GITHUB_RUN_ATTEMPT"])
        require(re.fullmatch(r"[0-9a-f]{40}", head) and run_id > 0 and attempt > 0,
                "invalid tag run identity")
        return cls(head, run_id, attempt)


def release_files(packet, receipt_digest, run):
    require({path.name for path in packet.iterdir()} == FILES | {"STAGED.json"},
            "staged packet must contain exactly six assets and the receipt")
    raw = (packet / "STAGED.json").read_bytes()
    require(sha256(raw) == receipt_digest, "stage receipt differs from the audit job output")
    receipt = json.loads(raw)
    context = {
        "schema": "rsomics.native-release-stage.v1", "repository": REPOSITORY,
        "product": PRODUCT, "version": VERSION, "head": run.head,
        "run_id": run.run_id, "attempt": run.attempt, "release_authorized": False,
    }
    require(all(receipt.get(key) == value for key, value in context.items()),
            "stage receipt is not from the exact tag run")
    files = {name: (packet / name).read_bytes() for name in sorted(FILES)}
    require(receipt["assets"] == {name: {"bytes": len(body), "sha256": sha256(body)}
                                  for name, body in files.items()}, "staged asset bytes differ")
    sums = "".join(f"{sha256(body)}  {name}\n" for name, body in files.items() if name != "SHA256SUMS")
    require(files["SHA256SUMS"] == sums.encode(), "staged checksums differ")
    return files


class GitHub:
    def __init__(self, token):
        require(token, "release job token is missing")
        self.token = token

    def __call__(self, method, path, payload=None, *, upload=False):
        host = "api.github.com"
        if upload:
            url = urlsplit(path)
            require(url.scheme == "https" and url.netloc == UPLOAD_HOST, "invalid asset upload host")
            host, path = UPLOAD_HOST, f"{url.path}?{url.query}"
        headers = {
            "Authorization": f"Bearer {self.token}", "Accept": "application/vnd.github+json",
            "X-GitHub-Api-Version": "2026-03-10", "User-Agent": f"{PRODUCT}-release",
        }
        body = None
        if payload is not None:
            headers["Content-Type"] = "application/octet-stream" if upload else "application/json"
            body = payload if upload else json.dumps(payload).encode()
        connection = HTTPSConnection(host, timeout=60)
        try:
            connection.request(method, path, body=body, headers=headers)
            response = connection.getresponse()
            expected_status = 201 if method == "POST" else 200
            require(response.status == expected_status,
                    f"GitHub {method} {path} returned HTTP {response.status}")
            return json.loads(response.read())
        finally:
            connection.close()


def pages(api, path):
    page = 1
    while True:
        rows = api("GET", f"{path}?per_page=100&page={page}")
        require(isinstance(rows, list), "expected a GitHub listing")
        yield from rows
        if len(rows) < 100:
            return
        page += 1


def verify_refs(api, head):
    for ref in ("heads/main", f"tags/{TAG}"):
        obj = api("GET", f"{API_ROOT}/git/ref/{ref}")["object"]
        require(obj["type"] == "commit" and obj["sha"] == head,
                f"remote {ref} does not point to the qualified commit")


def verify_release(record, release_id, notes, draft):
    require(record["id"] == release_id and record["tag_name"] == TAG
            and record["body"] == notes and record["draft"] is draft
            and record["prerelease"] is False, "GitHub release identity or state differs")


def verify_assets(records, files):
    observed, ids = {}, set()
    for asset in records:
        asset_id, name = asset["id"], asset["name"]
        require(type(asset_id) is int and asset_id > 0 and asset_id not in ids
                and name not in observed and asset["state"] == "uploaded",
                "duplicate or incomplete uploaded asset")
        ids.add(asset_id)
        observed[name] = asset["size"], asset["digest"]
    require(observed == {name: (len(body), f"sha256:{sha256(body)}") for name, body in files.items()},
            "uploaded asset names, sizes, or SHA-256 digests differ")


def publish(packet, receipt_digest, notes, env, event, api, report):
    run = TagRun.read(env, event)
    files = release_files(packet, receipt_digest, run)
    require(notes.strip(), "release notes are empty")
    verify_refs(api, run.head)
    require(not any(row["tag_name"] == TAG for row in pages(api, f"{API_ROOT}/releases")),
            f"a release or draft for {TAG} already exists")

    draft = api("POST", f"{API_ROOT}/releases", {
        "tag_name": TAG, "target_commitish": run.head, "name": f"{PRODUCT} {TAG}",
        "body": notes, "draft": True, "prerelease": False,
    })
    release_id = draft["id"]
    require(type(release_id) is int and release_id > 0, "invalid GitHub release ID")
    report(f"Draft release ID: {release_id}")
    verify_release(draft, release_id, notes, True)
    path = f"{API_ROOT}/releases/{release_id}"
    upload_url = f"https://{UPLOAD_HOST}{path}/assets"
    require(draft["upload_url"] in (upload_url, upload_url + "{?name,label}"),
            "unexpected GitHub asset upload URL")
    for name, body in files.items():
        api("POST", f"{upload_url}?name={quote(name, safe='')}", body, upload=True)

    verify_assets(pages(api, f"{path}/assets"), files)
    verify_release(api("GET", path), release_id, notes, True)
    verify_refs(api, run.head)
    released = api("PATCH", path, {"draft": False})
    verify_release(released, release_id, notes, False)
    verify_assets(released["assets"], files)
    return released


def report(message):
    print(message, flush=True)
    with Path(os.environ["GITHUB_STEP_SUMMARY"]).open("a") as summary:
        summary.write(message + "\n")


def main():
    parser = argparse.ArgumentParser(description="Publish the qualified FASTQ release from its new tag run")
    parser.add_argument("packet", type=Path)
    parser.add_argument("notes", type=Path)
    args = parser.parse_args()
    released = publish(args.packet, os.environ["RELEASE_STAGE_SHA256"], args.notes.read_text().strip(),
                        os.environ, json.loads(Path(os.environ["GITHUB_EVENT_PATH"]).read_bytes()),
                        GitHub(os.environ["RELEASE_TOKEN"]), report)
    report(f"Published release {released['id']}: https://github.com/{REPOSITORY}/releases/tag/{TAG}")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, KeyError, OSError, HTTPException) as error:
        message = str(error)[-2000:].replace("%", "%25").replace("\r", "%0D").replace("\n", "%0A")
        print(f"::error title=FASTQ release gate::{message}")
        raise SystemExit(1) from None
