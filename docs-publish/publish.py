"""Publish this run's finished HTML, then link the live site from the repository."""
import base64
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import time
import urllib.request

REPO = 'PRECISION-GM/BreezeLab.jl'
RUN = Path("/shared/home/greg/breezelab-work/integration")
BUILD = RUN / 'docs/build'

def command(args, **kwargs):
    return subprocess.run(args, check=True, text=True, capture_output=True, **kwargs).stdout.strip()

def api(path, method='GET', body=None):
    args = ['gh', 'api', path, '--method', method]
    if body is not None:
        args += ['--input', '-']
    output = command(args, input=json.dumps(body) if body is not None else None)
    return json.loads(output) if output else None

assert (BUILD / 'index.html').is_file()
assert (BUILD / 'generated/eastern_north_atlantic.png').is_file()
assert not list(BUILD.rglob('*.jld2')), 'Raw simulation data must not be published'
metadata = json.loads((RUN / 'gpu-build.json').read_text())
(BUILD / 'gpu-build.json').write_text(json.dumps(metadata, indent=2) + '\n')
(BUILD / '.nojekyll').touch()
identity = api('user')['login']
remote = f'https://github.com/{REPO}.git'
with tempfile.TemporaryDirectory(prefix='publish-', dir=RUN) as directory:
    directory = Path(directory)
    command(['git', 'init', '--initial-branch=gh-pages', str(directory)])
    def git(*args):
        return command(['git', '-C', str(directory), *args])
    git('remote', 'add', 'origin', remote)
    refs = git('ls-remote', '--heads', 'origin', 'gh-pages')
    if refs:
        git('fetch', 'origin', 'gh-pages')
        git('reset', '--hard', 'FETCH_HEAD')
        # Remove only this disposable publication checkout's old site files.
        for item in directory.iterdir():
            if item.name != '.git':
                shutil.rmtree(item) if item.is_dir() else item.unlink()
    shutil.copytree(BUILD, directory, dirs_exist_ok=True)
    git('config', 'user.name', identity)
    git('config', 'user.email', f'{identity}@users.noreply.github.com')
    git('add', '--all')
    git('commit', '-m', f"Publish GPU documentation from {metadata['source_commit'][:7]} (Slurm {metadata['job_id']})")
    pages_commit = git('rev-parse', 'HEAD')
    git('-c', 'credential.helper=', '-c', 'credential.helper=!gh auth git-credential', 'push', 'origin', 'gh-pages')

try:
    page = api(f'repos/{REPO}/pages')
except subprocess.CalledProcessError as error:
    if '404' not in error.stderr:
        raise
    page = api(f'repos/{REPO}/pages', 'POST',
               {'build_type': 'legacy', 'source': {'branch': 'gh-pages', 'path': '/'}})
url = page['html_url']
api(f'repos/{REPO}/pages/builds', 'POST')
for attempt in range(60):
    latest = api(f'repos/{REPO}/pages/builds/latest')
    if latest.get('commit') == pages_commit:
        if latest['status'] == 'errored':
            raise RuntimeError(latest)
        if latest['status'] == 'built':
            try:
                # Check the deployed build record, not merely an old index response.
                with urllib.request.urlopen(url + 'gpu-build.json', timeout=30) as response:
                    live = json.load(response)
                if live == metadata:
                    break
            except Exception:
                pass
    time.sleep(20)
else:
    raise RuntimeError('Pages has not served this GPU build after 20 minutes')

api(f'repos/{REPO}', 'PATCH', {'homepage': url})
readme_commit = None
for attempt in range(3):
    readme = api(f'repos/{REPO}/contents/README.md?ref=main')
    text = base64.b64decode(readme['content']).decode()
    if f']({url})' in text:
        break
    heading, rest = text.split('\n', 1)
    text = heading + '\n\n[Documentation](' + url + ')\n' + rest
    try:
        result = api(f'repos/{REPO}/contents/README.md', 'PUT',
                     {'message': 'Link GPU-generated documentation', 'branch': 'main',
                      'sha': readme['sha'], 'content': base64.b64encode(text.encode()).decode()})
        readme_commit = result['commit']['sha']
        break
    except subprocess.CalledProcessError as error:
        if attempt == 2 or '409' not in error.stderr:
            raise
else:
    raise RuntimeError('Could not update README')

result = dict(metadata, url=url, pages_commit=pages_commit, readme_commit=readme_commit)
(RUN / 'publication.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2), flush=True)
