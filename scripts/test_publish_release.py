import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from urllib.parse import parse_qs, urlsplit


REPO = "omics-rust/rsomics-fastq-preprocess"
ROOT = f"/repos/{REPO}"
REF = "refs/tags/v0.2.0"
HEAD = "1" * 40
NOTES = "FASTQ 0.2.0 release notes"
DRAFT = f"{ROOT}/releases/42"
UPLOAD = f"https://uploads.github.com{DRAFT}/assets"
FILES = {
    "rsomics-fastq-preprocess-0.2.0-linux-aarch64.tar.gz": b"Linux ARM bytes",
    "rsomics-fastq-preprocess-0.2.0-linux-x86_64.tar.gz": b"Linux Intel bytes",
    "rsomics-fastq-preprocess-0.2.0-macos-aarch64.tar.gz": b"macOS ARM bytes",
    "rsomics-fastq-preprocess-0.2.0-macos-x86_64.tar.gz": b"macOS Intel bytes",
    "rsomics-fastq-preprocess-0.2.0.crate": b"source package bytes",
}


def sha(body):
    return hashlib.sha256(body).hexdigest()


def make_packet(path):
    checksums = "".join(f"{sha(body)}  {name}\n" for name, body in sorted(FILES.items()))
    files = {**FILES, "SHA256SUMS": checksums.encode()}
    for name, body in files.items():
        (path / name).write_bytes(body)
    receipt = {
        "schema": "rsomics.native-release-stage.v1",
        "repository": REPO,
        "product": "rsomics-fastq-preprocess",
        "version": "0.2.0",
        "head": HEAD,
        "run_id": 123,
        "attempt": 1,
        "audit": {"bytes": 100, "sha256": "2" * 64},
        "release_authorized": False,
        "assets": {name: {"bytes": len(body), "sha256": sha(body)} for name, body in files.items()},
    }
    data = json.dumps(receipt).encode()
    (path / "STAGED.json").write_bytes(data)
    return files, sha(data)


class GitHubFixture:
    def __init__(self):
        self.calls = []
        self.uploaded = {}
        self.draft = None
        self.release_pages = [[]]
        self.tag_head = HEAD
        self.moved_after_upload = False
        self.bad_digest = False

    def assets(self):
        return [
            {"id": index, "name": name, "state": "uploaded", "size": len(body),
             "digest": "sha256:" + ("0" * 64 if self.bad_digest else sha(body))}
            for index, (name, body) in enumerate(sorted(self.uploaded.items()), start=100)
        ]

    def __call__(self, method, path, payload=None, *, upload=False):
        self.calls.append((method, path))
        if method == "GET" and path == f"{ROOT}/git/ref/heads/main":
            return {"object": {"type": "commit", "sha": HEAD}}
        if method == "GET" and path == f"{ROOT}/git/ref/tags/v0.2.0":
            head = "9" * 40 if self.moved_after_upload and self.uploaded else self.tag_head
            return {"object": {"type": "commit", "sha": head}}
        if method == "GET" and path.startswith(f"{ROOT}/releases?per_page=100&page="):
            return self.release_pages[int(path.rsplit("=", 1)[1]) - 1]
        if method == "POST" and path == f"{ROOT}/releases":
            assert payload == {"tag_name": "v0.2.0", "target_commitish": HEAD,
                               "name": "rsomics-fastq-preprocess v0.2.0", "body": NOTES,
                               "draft": True, "prerelease": False}
            self.draft = {**payload, "id": 42, "upload_url": UPLOAD + "{?name,label}"}
            return dict(self.draft)
        if method == "POST" and upload and path.startswith(UPLOAD + "?name="):
            assert self.draft["draft"]
            name, = parse_qs(urlsplit(path).query)["name"]
            assert name not in self.uploaded
            self.uploaded[name] = payload
            return next(asset for asset in self.assets() if asset["name"] == name)
        if method == "GET" and path == f"{DRAFT}/assets?per_page=100&page=1":
            return self.assets()
        if method == "GET" and path == DRAFT:
            return dict(self.draft)
        if method == "PATCH" and path == DRAFT:
            assert payload == {"draft": False}
            self.draft["draft"] = False
            return {**self.draft, "assets": self.assets(),
                    "html_url": f"https://github.com/{REPO}/releases/tag/v0.2.0"}
        raise AssertionError((method, path))


class PublishReleaseTests(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(importlib.util.find_spec("publish_release"), "release gate is missing")
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.packet = Path(directory.name)
        self.files, self.receipt_digest = make_packet(self.packet)
        self.env = {
            "GITHUB_EVENT_NAME": "push", "GITHUB_REF": REF, "GITHUB_REPOSITORY": REPO,
            "GITHUB_SHA": HEAD, "GITHUB_RUN_ID": "123", "GITHUB_RUN_ATTEMPT": "1",
        }
        self.event = {"ref": REF, "created": True, "deleted": False}
        self.api = GitHubFixture()
        self.reports = []

    def run_gate(self):
        from publish_release import publish

        return publish(self.packet, self.receipt_digest, NOTES, self.env, self.event,
                       self.api, self.reports.append)

    def test_exact_six_files_are_uploaded_and_verified_before_publishing(self):
        published = self.run_gate()
        self.assertEqual(self.api.uploaded, self.files)
        self.assertEqual(published["id"], 42)
        self.assertIs(published["draft"], False)
        self.assertEqual(self.api.calls[-1], ("PATCH", DRAFT))
        self.assertIn(("GET", f"{DRAFT}/assets?per_page=100&page=1"), self.api.calls[:-1])
        self.assertTrue(any("42" in line for line in self.reports))

    def test_wrong_ref_or_run_cannot_make_requests(self):
        self.env["GITHUB_REF"] = "refs/heads/main"
        with self.assertRaises(ValueError):
            self.run_gate()
        self.env["GITHUB_REF"] = REF
        self.env["GITHUB_RUN_ATTEMPT"] = "2"
        with self.assertRaises(ValueError):
            self.run_gate()
        self.assertEqual(self.api.calls, [])

    def test_tampered_receipt_or_asset_cannot_make_requests(self):
        good_digest = self.receipt_digest
        self.receipt_digest = "0" * 64
        with self.assertRaises(ValueError):
            self.run_gate()
        self.receipt_digest = good_digest
        (self.packet / "rsomics-fastq-preprocess-0.2.0.crate").write_bytes(b"different")
        with self.assertRaises(ValueError):
            self.run_gate()
        self.assertEqual(self.api.calls, [])

    def test_existing_draft_on_later_page_is_not_reused(self):
        self.api.release_pages = [[{"tag_name": f"old-{i}"} for i in range(100)],
                                  [{"tag_name": "v0.2.0", "id": 71, "draft": True}]]
        with self.assertRaises(ValueError):
            self.run_gate()
        self.assertTrue(all(method == "GET" for method, _ in self.api.calls))

    def test_wrong_initial_tag_head_cannot_create_draft(self):
        self.api.tag_head = "9" * 40
        with self.assertRaises(ValueError):
            self.run_gate()
        self.assertIsNone(self.api.draft)

    def test_wrong_digest_leaves_draft(self):
        self.api.bad_digest = True
        with self.assertRaises(ValueError):
            self.run_gate()
        self.assertEqual(self.api.uploaded, self.files)
        self.assertIs(self.api.draft["draft"], True)

    def test_moved_tag_after_upload_leaves_draft(self):
        self.api.moved_after_upload = True
        with self.assertRaises(ValueError):
            self.run_gate()
        self.assertEqual(self.api.uploaded, self.files)
        self.assertIs(self.api.draft["draft"], True)
