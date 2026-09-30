# Offline viewer assembly and local preview

Embeds archived presentation arrays and repository-native scripts into a single
HTML artifact. `build` does not modify the simulation or its results. The optional
preview server binds only localhost and serves only the showcase output folder.
Its one write endpoint accepts a bounded PNG payload at a fixed poster path;
it accepts no arbitrary filesystem paths. Stop the server when preview work
ends; the standalone HTML remains usable offline.
