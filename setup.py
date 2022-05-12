
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
    def compile_seq():
      target_dir = os.path.join(self.build_lib, 'biser/exe')
      self.mkpath(target_dir)
      env = os.environ.copy()
      if shutil.which("seqc"):
        env["SEQ_LIBRARY_PATH"] = Path(os.path.dirname(shutil.which("seqc"))) / ".." / "lib" / "seq"
      subprocess.check_call(
        ['seqc', 'build', 'biser/seq/__init__.seq', '-release', '-o', f'{target_dir}/biser.exe'],
        env=env
      )
      if sys.platform == "darwin":
        ext = "dylib"
        subprocess.check_call(
          ["install_name_tool", "-add_rpath", "@executable_path/.", f'{target_dir}/biser.exe']
        )
      else:
        ext = "so"
        subprocess.check_call(
          ["patchelf", "--set-rpath", "$ORIGIN", f'{target_dir}/biser.exe']
        )
      os.unlink(f'{target_dir}/biser.exe.o')
      seqpath = Path(shutil.which('seqc')).parent
      for lib in ["libseqrt", "libomp"]:
        for p in [
          seqpath / f"{lib}.{ext}",
          seqpath / ".." / "lib" / "seq" / f"{lib}.{ext}"
        ]:
          print(p)
          if os.path.exists(p):
            subprocess.check_call(["cp", p, f'{target_dir}/{lib}.{ext}'])
            break
    if not self.dry_run:
      compile_seq()
    build_py.run(self)


exec(open("biser/version.py").read())
setup(
  name="biser",
  python_requires='>=3.7',
  version=__version__,
  description="Fast Characterization of Segmental Duplication Structure in Multiple Genome Assemblies",
  url="https://github.com/0xTCG/biser",
  long_description="Please see https://github.com/0xTCG/biser for more details.",
  author="Hamza Išerić, Ibrahim Numanagić",
  author_email="inumanag@uvic.ca",
  download_url="https://github.com/0xTCG/biser/tarball/master",
  license="MIT License.",
  keywords=["genome analysis", "fast alignment", "segmental duplications", "sequence decomposition"],
  install_requires=["tqdm", "ncls"],
  entry_points={"console_scripts": ["biser = biser.__main__:console"]},
  packages=find_packages(),
  package_data={ "biser": ["biser/seq"], },
  include_package_data=True,
  cmdclass={
    'build_py': CustomBuild,
  },
)