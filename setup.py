from setuptools import setup

setup(
   name='DynCode',
   version='1',
   description='Time series analysis package',
   author='fiorefabris',
   author_email='fiorefabris@gmail.com.com',
   packages=['DynCode', 'DynCode.LocalAnalysis','DynCode.DataAnalysis'],
   install_requires=['matplotlib', 'pandas','numpy'], #external packages as dependencies
)

Python 3.4 or higher!
