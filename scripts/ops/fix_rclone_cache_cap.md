# The OneDrive mount has no VFS cache cap — one line fixes it

**Status: NOT APPLIED. Needs a remount, which would interrupt another project's in-flight
writes. Flagged for a moment when the box is quiet.**

## What was found, 2026-09-17

The watchdog reported the rclone VFS cache at **5.51 GB**, up from 0.0 GB on every previous
tick. Verified by two methods because they disagreed:

| method | says |
|---|---|
| `du -sh C:/Temp/rclone-cache` | 192 MB |
| Python `stat()` over 17 files | **5.51 GB**, block count matching apparent size |

`du` under Git Bash is the one that is wrong here — it under-reports these files on NTFS.
The watchdog is correct. (My first reading was that the watchdog had a sparse-file bug; it
does not. Measuring twice is what caught it.)

## What is in the cache

| size | file | date |
|---|---|---|
| 5.06 GB | `OpenDEL_libraries.zip` | 2026-07-11 |
| ~0.03 GB × 16 | `cognitive-hub-v*.apk` | rolling, latest 2026-09-17 00:12 |

**None of it belongs to this campaign.** The APKs are the edgemere-mobile project; the zip
is a DEL library archive. Not ours to delete, and the ops rule is explicit about that.

## Why it matters anyway

`C:\Temp\tools\mount-onedrive.vbs` mounts with `--vfs-cache-mode full` and **no
`--vfs-cache-max-size`**. Every byte read *or written* through `O:\` is cached on C:
indefinitely and never evicted. This repo's CLAUDE.md records it growing to 8 GB once and
driving C: to zero, which broke unrelated tools on the box.

C: currently has ~10 GB free, so this is not urgent. It is also not self-limiting: nothing
in the current configuration will ever shrink that cache.

## The fix

Add two flags to the mount command in `C:\Temp\tools\mount-onedrive.vbs`:

```
--vfs-cache-max-size 8G --vfs-cache-max-age 24h
```

Then remount (`log off/on`, or kill the rclone process and re-run the .vbs).

## Why it was not applied here

A remount drops the mount briefly. Something was writing cognitive-hub APKs through it as
recently as 00:12 today, and interrupting another agent's upload to save disk that is not
currently scarce is a bad trade. This wants a quiet moment and a deliberate decision, not
an autonomous ops tick.

**Not a workaround:** this campaign already avoids the hazard entirely — `cypstruct.storage`
pushes through the `onedrive:` *remote* via the Graph API, which never touches the VFS
cache. That is why our own 1.7 GiB archive left the cache at 0.0 GB. The fix is for
everything else on the box.
