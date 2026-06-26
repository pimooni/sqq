from copy import deepcopy
import unittest

from sqq.cli import build_parser
from sqq.config import DEFAULT_CONFIG
from sqq.pipeline import normalize_analysis_scopes, resolve_cage_report_types


def normalized_config(search_sizes, report_types="auto"):
    config = deepcopy(DEFAULT_CONFIG)
    config["ring"]["sizes"] = list(search_sizes)
    config["cage"]["report_types"] = report_types
    normalize_analysis_scopes(config)
    return config


class CageReportScopeTests(unittest.TestCase):
    def test_cli_help_uses_public_cage_groups(self):
        parser = build_parser()
        subparsers = next(action for action in parser._actions if getattr(action, "dest", None) == "command")
        help_text = subparsers.choices["analyze"].format_help()

        for value in ("auto", "all", "I", "II", "H", "HS-I", "TS-I", "I2II"):
            self.assertIn(value, help_text)
        self.assertNotIn("51268", help_text)
        self.assertNotIn("435663", help_text)

    def test_default_cage_reports_follow_search_sizes(self):
        config = normalized_config([4, 5, 6])

        self.assertEqual(config["ring"]["sizes"], [4, 5, 6])
        self.assertEqual(config["cage"]["report_types"], "all")
        self.assertIsNone(resolve_cage_report_types("auto", [4, 5, 6], 20))

    def test_structure_groups_expand_and_deduplicate(self):
        self.assertEqual(
            resolve_cage_report_types("I,II", [4, 5, 6], 20),
            ("512", "51262", "51264"),
        )

    def test_h_group_expands_to_type_h_cages(self):
        self.assertEqual(
            resolve_cage_report_types("H", [4, 5, 6], 20),
            ("512", "51268", "435663"),
        )

    def test_transition_group_expands_to_51263(self):
        self.assertEqual(resolve_cage_report_types("I2II", [5, 6], 20), ("51263",))

    def test_exact_cage_report_types_remain_supported(self):
        config = normalized_config([4, 5, 6], ["512", "435663"])

        self.assertEqual(config["cage"]["report_types"], ["512", "435663"])

    def test_all_is_an_explicit_alias_for_auto_reporting(self):
        self.assertIsNone(resolve_cage_report_types("all", [4, 5, 6], 20))

    def test_group_members_must_fit_search_sizes(self):
        with self.assertRaisesRegex(ValueError, "absent from --size"):
            normalized_config([5, 6], "H")


if __name__ == "__main__":
    unittest.main()
