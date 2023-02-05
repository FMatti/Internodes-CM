# Internodes-CM

Implementation of the internodes method for two- and three-dimensional problems in contact mechanics.

## Info

This code served as a prototype for the inclusion of the method in [Akantu](https://gitlab.com/akantu/akantu/-/tree/features/40-contact-using-the-internodes-method) in collaboration with [@Technici4n](https://github.com/Technici4n).

## Instructions

Clone the repository
```[bash]
git clone https://github.com/FMatti/Internodes-CM.git
cd Internodes-CM
```

Optionally create and activate a virtual environment
```[bash]
python -m venv .venv

source .venv/bin/activate (on Linux, macOS)
.venv\Scripts\activate.bat (on Windows)
```

Install the required packages
```[bash]
pip install --upgrade pip
pip install -r requirements.txt
```

Reproduce our results (takes ~5 min)
```[bash]
python main.py
```

## 