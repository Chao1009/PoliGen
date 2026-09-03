Snippet for README.md's Build section (not applied there — that file is
owned by another running job; paste by hand when free):

```markdown
Alternatively, `pip install` the Python package on its own (scikit-build-core;
no `env.sh`, no `build/` needed):

    LIPOLGEN_DEPS_PREFIX=/path/to/deps/install pip install -e .

`PYTHIA8DATA`/`LHAPDF_DATA_PATH` still need exporting at run time (see
`docs/USAGE.md`); the wheel's RPATH points at this machine's deps prefix, so
it is not relocatable as-is — `auditwheel repair` fixes that but pulls in
GPL-3 redistribution terms for the combined work. See `docs/USAGE.md` for
details.
```
