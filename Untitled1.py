from setuptools import setup, find_packages

with open('old_requirements.txt') as f:
    requirements = f.read().splitlines()

for i in range(len(requirements)):
    requirements[i] = requirements[i].split(' ')

requirements = [item for sublist in requirements for item in sublist]
requirements = sorted(list(set(requirements)))

with open('requirements.txt', 'w') as fileOut:
    for i in range(len(requirements)):
        fileOut.write(requirements[i] + '\n')
