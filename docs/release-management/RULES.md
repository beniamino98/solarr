# Release Management Rules

## Versioning

Production versions use three numeric components: `MAJOR.MINOR.PATCH`.

- `PATCH`: compatible bug fix, documentation update, test update, or low-risk UI correction. Example: `0.0.0` -> `0.0.1`.
- `MINOR`: compatible feature, workflow, validation, import, or export. Example: `0.0.1` -> `0.1.0`.
- `MAJOR`: incompatible change requested or accepted by the user. Example: `0.4.2` -> `1.0.0`.

## Branch Workflow

1. Add or edit exactly one active issue in `TRACKER.md` before implementation.
2. Set the issue `Change type` to `PATCH`, `MINOR`, or `MAJOR`.
3. Saving the issue marks it as the only `In progress` issue.
4. Resolve the working branch from the issue change type and any existing release branch already present in `VERSIONING.md` or local Git.
5. Record the established working branch in `README.md` metadata.
6. Create the release branch only when no suitable existing working branch exists.
7. Use branch names in this format:

```text
PREFIX-TARGET_VERSION
```

Example: `PB-0.4.2`.

## Commit Workflow

Commit messages begin with the working branch target version and a sequential three-digit commit index.
The index is the next available number after existing commits already made in that branch.

```text
PREFIX-TARGET_VERSION-005 - ISSUE_ID
PREFIX-TARGET_VERSION-005 - ISSUE_ID; ISSUE_ID
```

Examples:

```text
PB-0.4.2-005 - PB-0020
PB-0.4.2-005 - PB-0020; PB-0032
```

## Push And Merge Workflow

1. Commit completed changes on the release branch.
2. Push the release branch when remote synchronization is required.
3. Use a structured production merge message in this form:

```text
PREFIX-TARGET_VERSION-MERGE - WORKING_BRANCH -> main - ISSUE_ID
```

4. Merge the release branch into `main` only after relevant issues are verified and closed into `VERSIONING.md`.
5. After merge, the target production version becomes the new current production version.

## General Rules

- Closed and verified issues are removed from `TRACKER.md` and recorded in `VERSIONING.md`.
- Only one open issue should be `In progress` at a time.
- The version where an issue was opened stays attached to the issue even if it is fixed in a later release.
- Branch/version warnings should be resolved before commit or merge.

## Custom Project Rules

_No custom project rules declared._
