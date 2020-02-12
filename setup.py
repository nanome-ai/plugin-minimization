import pathlib
from setuptools import find_packages, setup

README = (pathlib.Path(__file__).parent / "README.md").read_text()

setup(
    name = 'nanome-minimization',
    packages=find_packages(),
    version = '0.2.1',
    license='MIT',
    description = 'Nanome Plugin to run minimization on selected structures',
    long_description = README,
    long_description_content_type = "text/markdown",
    author = 'Nanome',
    author_email = 'hello@nanome.ai',
    url = 'https://github.com/nanome-ai/plugin-minimization',
    platforms="any",
    keywords = ['virtual-reality', 'chemistry', 'python', 'api', 'plugin'],
    install_requires=['nanome'],
    entry_points={"console_scripts": ["nanome-minimization = nanome_minimization.Minimization:main"]},
    classifiers=[
        'Development Status :: 3 - Alpha',

        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry',

        'License :: OSI Approved :: MIT License',

        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    package_data={
        "nanome_minimization": [
            "_MinimizationMenu.json"
        ]
    },
)
