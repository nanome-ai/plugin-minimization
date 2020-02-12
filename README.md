# Nanome - Minimization

Minimization arranges the selected molecule to try to find a geometry where inter-atomic forces are as close to zero as possible.

### Dependencies

`nanome-minimization` requires `nanobabel` to exist in the working directory, or for the path to be specified in the environment variable `NANOBABEL`.

### Installation

```sh
$ pip install nanome-minimization
```

### Usage

To start the plugin:

```sh
$ nanome-minimization -a plugin_server_address
```

In Nanome:

- Activate Plugin
- Select atoms to minimize
- Click Run to immediately minimize or Advanced Settings to choose different options

### License

MIT
