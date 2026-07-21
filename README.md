# BoolDog <img src="docs/figures/logo.png" raw=true alt="BoolDoG icon"  width="240" align="right" >

[![PyPI](https://img.shields.io/pypi/v/booldog)](https://pypi.org/project/booldog/)
[![Python >=3.12](https://img.shields.io/badge/python-%3E%3D3.12-blue)](https://pypi.org/project/booldog/)
[![License: GPL-3.0](https://img.shields.io/github/license/NIB-SI/BoolDog)](LICENSE)
[![Docs](https://img.shields.io/github/actions/workflow/status/NIB-SI/BoolDog/deploy-docs.yml?label=docs)](https://nib-si.github.io/BoolDog)
[![Tests](https://img.shields.io/github/actions/workflow/status/NIB-SI/BoolDog/tests.yml?label=tests)](https://github.com/NIB-SI/BoolDog/actions/workflows/tests.yml)

A Python package for integrated Boolean and semi-quantitative network modelling.

## CI Status

**Tests**

| OS \ Python | 3.12 | 3.13 |
|---|---|---|
| Ubuntu | [![Ubuntu 3.12](https://img.shields.io/github/check-runs/NIB-SI/BoolDog/main?nameFilter=test%20%28ubuntu-latest%2C%203.12%29)](https://github.com/NIB-SI/BoolDog/actions/workflows/tests.yml) | [![Ubuntu 3.13](https://img.shields.io/github/check-runs/NIB-SI/BoolDog/main?nameFilter=test%20%28ubuntu-latest%2C%203.13%29)](https://github.com/NIB-SI/BoolDog/actions/workflows/tests.yml) |
| macOS | [![macOS 3.12](https://img.shields.io/github/check-runs/NIB-SI/BoolDog/main?nameFilter=test%20%28macos-latest%2C%203.12%29)](https://github.com/NIB-SI/BoolDog/actions/workflows/tests.yml) | [![macOS 3.13](https://img.shields.io/github/check-runs/NIB-SI/BoolDog/main?nameFilter=test%20%28macos-latest%2C%203.13%29)](https://github.com/NIB-SI/BoolDog/actions/workflows/tests.yml) |
| Windows | [![Windows 3.12](https://img.shields.io/github/check-runs/NIB-SI/BoolDog/main?nameFilter=test%20%28windows-latest%2C%203.12%29)](https://github.com/NIB-SI/BoolDog/actions/workflows/tests.yml) | [![Windows 3.13](https://img.shields.io/github/check-runs/NIB-SI/BoolDog/main?nameFilter=test%20%28windows-latest%2C%203.13%29)](https://github.com/NIB-SI/BoolDog/actions/workflows/tests.yml) |

**Tutorials**

<table>
<tr><td valign="top">
<b>tutorial-basic</b><br>
<table>
<tr><th>OS \ Python</th><th>3.12</th><th>3.13</th></tr>
<tr><td>Ubuntu</td><td><a href="https://github.com/NIB-SI/BoolDog/actions/workflows/tests.yml"><img src="https://img.shields.io/github/check-runs/NIB-SI/BoolDog/main?nameFilter=tutorials%20%28ubuntu-latest%2C%203.12%2C%20tutorial-basic%29" alt="Ubuntu 3.12"></a></td><td><a href="https://github.com/NIB-SI/BoolDog/actions/workflows/tests.yml"><img src="https://img.shields.io/github/check-runs/NIB-SI/BoolDog/main?nameFilter=tutorials%20%28ubuntu-latest%2C%203.13%2C%20tutorial-basic%29" alt="Ubuntu 3.13"></a></td></tr>
<tr><td>macOS</td><td><a href="https://github.com/NIB-SI/BoolDog/actions/workflows/tests.yml"><img src="https://img.shields.io/github/check-runs/NIB-SI/BoolDog/main?nameFilter=tutorials%20%28macos-latest%2C%203.12%2C%20tutorial-basic%29" alt="macOS 3.12"></a></td><td><a href="https://github.com/NIB-SI/BoolDog/actions/workflows/tests.yml"><img src="https://img.shields.io/github/check-runs/NIB-SI/BoolDog/main?nameFilter=tutorials%20%28macos-latest%2C%203.13%2C%20tutorial-basic%29" alt="macOS 3.13"></a></td></tr>
<tr><td>Windows</td><td><a href="https://github.com/NIB-SI/BoolDog/actions/workflows/tests.yml"><img src="https://img.shields.io/github/check-runs/NIB-SI/BoolDog/main?nameFilter=tutorials%20%28windows-latest%2C%203.12%2C%20tutorial-basic%29" alt="Windows 3.12"></a></td><td><a href="https://github.com/NIB-SI/BoolDog/actions/workflows/tests.yml"><img src="https://img.shields.io/github/check-runs/NIB-SI/BoolDog/main?nameFilter=tutorials%20%28windows-latest%2C%203.13%2C%20tutorial-basic%29" alt="Windows 3.13"></a></td></tr>
</table>
</td><td valign="top">
<b>tutorial-advanced</b><br>
<table>
<tr><th>OS \ Python</th><th>3.12</th><th>3.13</th></tr>
<tr><td>Ubuntu</td><td><a href="https://github.com/NIB-SI/BoolDog/actions/workflows/tests.yml"><img src="https://img.shields.io/github/check-runs/NIB-SI/BoolDog/main?nameFilter=tutorials%20%28ubuntu-latest%2C%203.12%2C%20tutorial-advanced%29" alt="Ubuntu 3.12"></a></td><td><a href="https://github.com/NIB-SI/BoolDog/actions/workflows/tests.yml"><img src="https://img.shields.io/github/check-runs/NIB-SI/BoolDog/main?nameFilter=tutorials%20%28ubuntu-latest%2C%203.13%2C%20tutorial-advanced%29" alt="Ubuntu 3.13"></a></td></tr>
<tr><td>macOS</td><td><a href="https://github.com/NIB-SI/BoolDog/actions/workflows/tests.yml"><img src="https://img.shields.io/github/check-runs/NIB-SI/BoolDog/main?nameFilter=tutorials%20%28macos-latest%2C%203.12%2C%20tutorial-advanced%29" alt="macOS 3.12"></a></td><td><a href="https://github.com/NIB-SI/BoolDog/actions/workflows/tests.yml"><img src="https://img.shields.io/github/check-runs/NIB-SI/BoolDog/main?nameFilter=tutorials%20%28macos-latest%2C%203.13%2C%20tutorial-advanced%29" alt="macOS 3.13"></a></td></tr>
<tr><td>Windows</td><td><a href="https://github.com/NIB-SI/BoolDog/actions/workflows/tests.yml"><img src="https://img.shields.io/github/check-runs/NIB-SI/BoolDog/main?nameFilter=tutorials%20%28windows-latest%2C%203.12%2C%20tutorial-advanced%29" alt="Windows 3.12"></a></td><td><a href="https://github.com/NIB-SI/BoolDog/actions/workflows/tests.yml"><img src="https://img.shields.io/github/check-runs/NIB-SI/BoolDog/main?nameFilter=tutorials%20%28windows-latest%2C%203.13%2C%20tutorial-advanced%29" alt="Windows 3.13"></a></td></tr>
</table>
</td></tr>
</table>

## Documentation

For installation, usage, and examples, see the documentation at 
[nib-si.github.io/BoolDog](https://nib-si.github.io/BoolDog)

## Citation

Bleker, C., Zagorščak, M., Blejec, A., Gruden, K. & Županič, A. BoolDog: integrated
Boolean and semi-quantitative network modelling in Python. *bioRxiv*
[10.64898/2026.03.16.711264](https://www.biorxiv.org/content/10.64898/2026.03.16.711264v1) (2026).

```bibtex
@article{bleker2026booldog,
  title   = {BoolDog: integrated Boolean and semi-quantitative network modelling in Python},
  author  = {Bleker, Carissa and Zagor{\v{s}}{\v{c}}ak, Maja and Blejec, Andrej and Gruden, Kristina and {\v{Z}}upani{\v{c}}, An{\v{z}}e},
  journal = {bioRxiv},
  year    = {2026},
  doi     = {10.64898/2026.03.16.711264}
}
```

## Development version

Development version of BoolDog can be installed from GitHub.

To install:

```bash
git clone https://github.com/NIB-SI/BoolDog.git
cd BoolDoG
pip install .
```

To install with all optional extras (networks, SBML-qual, BioModels):

```bash
pip install .[all]
```

To remove:

```bash
pip uninstall booldog
```

## Docker

A standalone [`Dockerfile`](Dockerfile) at the repo root builds a container
with BoolDog and its `all` extras installed. For a container that instead
extends the [CoLoMoTo Docker](https://github.com/colomoto/colomoto-docker)
image (bundled with its own tools and a Jupyter notebook environment), see
[colomoto/README.md](colomoto/README.md).

## Benchmarks

Performance/scalability characterisation of BoolDog runtime against model
size and regulatory in-degree is in
[benchmarks/README.md](benchmarks/README.md).

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md).

---

