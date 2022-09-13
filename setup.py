from setuptools import setup, find_packages
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    author="Jiaying Lai",
    description="A package for computing variation numbers",
    name="variation_number",
    version="0.2.1",
    packages=find_packages(include=["variation_number","variation_number*"]),
    install_requires=[
        'bio>=0.4.1',
        'numpy>=1.18.5',
        'scikit-bio>=0.5.5',
        'beautifulsoup4>=4.7.1'
    ],
    python_requires='>=3.7',
    long_description=long_description,
    long_description_content_type='text/markdown'
)