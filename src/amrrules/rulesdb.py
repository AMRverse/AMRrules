import json
import os
import shutil
import tarfile
import tempfile
import urllib.error
import urllib.request
from pathlib import Path

RULES_DB_REPO = "AMRverse/AMRrulesDB"
RULES_DB_ASSET = "amrrulesdb.tar.gz"
# bump this when the engine can no longer read rule files written to an older schema
DB_SCHEMA_VERSION = 1


def get_data_dir():
    """Where downloaded resources live. Defaults to the installed package folder; override with AMRRULES_DATA_DIR."""
    env = os.environ.get("AMRRULES_DATA_DIR")
    return Path(env).expanduser() if env else Path(__file__).resolve().parent


def get_rules_dir():
    """Folder holding the installed rules database. Exits with a helpful message if it's missing."""
    rules_dir = get_data_dir() / "rules"
    if not (rules_dir / "rule_key_file.tsv").exists():
        raise SystemExit(
            "The AMRrules rules database has not been downloaded yet.\n"
            "Please run: amrrules --download-resources rules\n"
            "Then rerun your original command."
        )
    return rules_dir


def get_rules_db_version():
    manifest = get_data_dir() / "rules" / "manifest.json"
    if not manifest.exists():
        return "unknown"
    return json.loads(manifest.read_text()).get("db_version", "unknown")


def rules_db_url(version=None):
    base = f"https://github.com/{RULES_DB_REPO}/releases"
    if version:
        return f"{base}/download/{version}/{RULES_DB_ASSET}"
    return f"{base}/latest/download/{RULES_DB_ASSET}"


def download_rules_db(version=None, url=None):
    """Download a rules database release, check it matches this engine, and swap it in."""
    url = url or rules_db_url(version)
    data_dir = get_data_dir()
    data_dir.mkdir(parents=True, exist_ok=True)

    print(f"Downloading rules database from {url}")
    # the temp folder lives inside data_dir so the final move is a same-filesystem rename
    with tempfile.TemporaryDirectory(dir=data_dir) as tmp:
        archive = Path(tmp) / RULES_DB_ASSET
        try:
            urllib.request.urlretrieve(url, archive)
        except urllib.error.URLError as exc:
            raise SystemExit(
                f"Could not download the rules database: {exc}\n"
                "Check your internet connection, and that the requested --db-version exists."
            ) from None

        staging = Path(tmp) / "rules"
        with tarfile.open(archive) as tar:
            tar.extractall(staging, filter="data")

        manifest = json.loads((staging / "manifest.json").read_text())
        if manifest.get("schema_version") != DB_SCHEMA_VERSION:
            raise SystemExit(
                f"Rules database {manifest.get('db_version')} uses schema {manifest.get('schema_version')}, "
                f"but this version of AMRrules expects schema {DB_SCHEMA_VERSION}. "
                "Update AMRrules, or choose a different --db-version."
            )

        rules_dir = data_dir / "rules"
        if rules_dir.exists():
            shutil.rmtree(rules_dir)
        shutil.move(str(staging), str(rules_dir))

    print(f"Rules database {manifest.get('db_version')} installed at {rules_dir}")