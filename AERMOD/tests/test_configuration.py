"""Check Python syntax, configuration bindings and notebook cleanliness."""
from pathlib import Path
import ast,json,unittest

ROOT=Path(__file__).resolve().parents[1]


class ConfigurationTest(unittest.TestCase):
    def test_script_configuration(self):
        cfg=json.loads((ROOT/'config/settings.example.json').read_text())
        for file in ROOT.rglob('*.py'):
            code=file.read_text(encoding='utf-8');tree=ast.parse(code,feature_version=(3,11))
            for node in ast.walk(tree):
                if isinstance(node,ast.Call) and isinstance(node.func,ast.Name) and node.args and isinstance(node.args[0],ast.Constant):
                    name=node.func.id;key=node.args[0].value
                    if name in {'_configured_path','path'}:self.assertIn(key,cfg['paths'],str(file))
                    if name=='_setting':self.assertIn(key,cfg['settings'],str(file))
        self.assertTrue(all(not Path(v).is_absolute() for v in cfg['paths'].values()))
    def test_notebook(self):
        notebook=json.loads((ROOT/'notebooks/Report_figures.ipynb').read_text(encoding='utf-8'))
        for cell in notebook['cells']:
            if cell['cell_type']=='code':
                tree=ast.parse(''.join(cell['source']),feature_version=(3,11))
                for node in tree.body:
                    if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='original_plot_code' for t in node.targets):
                        ast.parse(ast.literal_eval(node.value),feature_version=(3,11))
                self.assertEqual(cell['outputs'],[])
                self.assertIsNone(cell['execution_count'])


if __name__=='__main__':unittest.main()
