"""Reject altered measurement dependencies and apply the child memory budget."""

import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from run_cli import digest, run_one, verify_lock


class SoftwareLockTests(unittest.TestCase):
    def test_binary_parameter_and_source_changes_are_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            cache, binary, parameter = (root / name for name in ("CMakeCache.txt", "program", "energy.par"))
            for path in (cache, binary, parameter):
                path.write_text(path.name)
            lock = {"schema_version": 1, "source_commit": "fixed",
                    "cmake_cache_sha256": digest(cache),
                    "energy_parameter_sha256": digest(parameter),
                    "binaries": {"program": digest(binary)}}
            lock_path = root / "lock.json"
            lock_path.write_text(json.dumps(lock))
            with patch("run_cli.git_state", return_value={"commit": "fixed", "dirty": False}):
                self.assertEqual(verify_lock(lock_path, root, root, [{"name": "program"}], parameter), lock)
                for path in (binary, parameter, cache):
                    original = path.read_bytes()
                    path.write_bytes(original + b"changed")
                    with self.assertRaises(ValueError):
                        verify_lock(lock_path, root, root, [{"name": "program"}], parameter)
                    path.write_bytes(original)
            with patch("run_cli.git_state", return_value={"commit": "fixed", "dirty": True}):
                with self.assertRaises(ValueError):
                    verify_lock(lock_path, root, root, [{"name": "program"}], parameter)

    def test_child_receives_address_space_limit(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            command = [sys.executable, "-c", "import resource; print(resource.getrlimit(resource.RLIMIT_AS))"]
            result = run_one(command, root, os.environ.copy(), 5, memory_limit_mib=128)
            self.assertEqual(result["status"], "success")
            limit = 128 * 1024 * 1024
            self.assertEqual((root / "stdout.txt").read_text().strip(), str((limit, limit)))


if __name__ == "__main__":
    unittest.main()
