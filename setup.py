#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
import os
import sys
from pathlib import Path
import shutil
import subprocess


class CustomBuild(build_py):
    def run(self):
        def compile_codon():
            target_dir = os.path.join(self.build_lib, "biser/exe")
            self.mkpath(target_dir)
            env = os.environ.copy()
            subprocess.check_call(
                [
                    "codon",
                    "build",
                    "-plugin",
                    "seq",
                    "biser/codon/__init__.codon",
                    "-release",
                    "-o",
                    f"{target_dir}/biser.exe",
                ],
                env=env,
            )
            if sys.platform == "darwin":
                ext = "dylib"
                subprocess.check_call(
                    [
                        "install_name_tool",
                        "-add_rpath",
                        "@executable_path/.",
                        f"{target_dir}/biser.exe",
                    ]
                )
            else:
                ext = "so"
                subprocess.check_call(
                    ["patchelf", "--set-rpath", "$ORIGIN", f"{target_dir}/biser.exe"]
                )
            codon_path = Path(shutil.which("codon")).parent
            for lib in ["libcodonrt", "libomp"]:
                for p in [
                    codon_path / f"{lib}.{ext}",
                    codon_path / ".." / "lib" / "codon" / f"{lib}.{ext}",
                ]:
                    print(p)
                    if os.path.exists(p):
                        subprocess.check_call(["cp", p, f"{target_dir}/{lib}.{ext}"])
                        break

        if not self.dry_run:
            compile_codon()
        build_py.run(self)


exec(open("biser/version.py").read())
setup(
    name="biser",
    python_requires=">=3.7",
    version=__version__,
    description="Fast Characterization of Segmental Duplication Structure in Multiple Genome Assemblies",
    url="https://github.com/0xTCG/biser",
    long_description="Please see https://github.com/0xTCG/biser for more details.",
    author="Hamza Išerić, Ibrahim Numanagić",
    author_email="inumanag@uvic.ca",
    download_url="https://github.com/0xTCG/biser/tarball/master",
    license="MIT License.",
    keywords=[
        "genome analysis",
        "fast alignment",
        "segmental duplications",
        "sequence decomposition",
    ],
    install_requires=["tqdm", "ncls", "multiprocess"],
    entry_points={"console_scripts": ["biser = biser.__main__:console"]},
    packages=find_packages(),
    cmdclass={"build_py": CustomBuild},
)
