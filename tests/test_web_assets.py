import pathlib
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]


class TestWebApplicationAssets(unittest.TestCase):
    def test_index_links_versioned_runtime_and_local_assets(self):
        html = (ROOT / "index.html").read_text(encoding="utf-8")
        self.assertIn("pyodide/v0.29.5/full/pyodide.js", html)
        self.assertIn("./assets/app.js", html)
        self.assertIn("./assets/styles.css", html)
        self.assertIn("Content-Security-Policy", html)
        self.assertIn("Exploratory use only", html)

    def test_browser_uses_python_engine_without_unsafe_dom_html_injection(self):
        script = (ROOT / "assets" / "app.js").read_text(encoding="utf-8")
        self.assertIn("./admet_predictor/__init__.py", script)
        self.assertIn("ADMETPredictor.evaluate_candidate", script)
        self.assertIn("PharmacokineticSimulator.simulate_oral_multiple", script)
        self.assertNotIn(".innerHTML", script)

    def test_pages_workflow_pins_action_commits(self):
        workflow = (ROOT / ".github" / "workflows" / "pages.yml").read_text(encoding="utf-8")
        self.assertIn("actions/deploy-pages@368f82528645a54fb793d4d04e342629a3f51346", workflow)
        self.assertNotIn("uses: actions/deploy-pages@v", workflow)


if __name__ == "__main__":
    unittest.main()
