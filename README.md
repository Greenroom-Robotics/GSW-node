# GSW-node

Nodejs interface to Gibbs-SeaWater (GSW) Oceanographic Toolbox in C++.

Link to the homepage: [Teos-10.org](https://www.teos-10.org/)

This library includes all modules:

- TeosBase
- TeosIce
- TeosSea

## Requirements

- Node: 20
- cmake
- gcc
- python 3

Currently, this has only been tested on linux Ubuntu 24.04.

## TeosNode

Install with

```bash
$ pnpm add @greenroom-robotics/gsw_cpp
$ pnpm add @greenroom-robotics/gsw_node
```

If you clone this repo you can install using:

```bash 
$ pnpm i
```

Then test with:

```
$ pnpm test
```

### Use

Avaliable functions are generated and avaliable in: [lib-types.d.ts](lib-types.d.ts)

The use of `gsw_z_from_p` is shown in the below example, this is an example
from https://www.teos-10.org/pubs/gsw/html/gsw_z_from_p.html

**Note** if you see the definition of `gsw_z_from_p` in the source code you will notice that there are two optional
arguments, in the bindings of this function these defaults do not transfer over, in the following typescript code you
can see that these defaults need to be added .

```c++
double gsw_z_from_p(double p, double lat, double geo_strf_dyn_height=0.0,
								  double sea_surface_geopotential=0.0);
```

```typescript
import {TeosBase} from '@greenroom-robotics/gsw_node'

const teosBase = new TeosBase()
const p = 10
const lat = 4
const expectedResult = -0.099445834469453 * 1.0e+002
const result = teosBase.gsw_z_from_p(p, lat, 0, 0)

```

## TeosCpp

This is a copy of the source code from Teos-10, downloadable from https://www.teos-10.org/software.htm

## Other notes

Not all functions have been implemented for node js.
To add new functions update [gsw-node](packages/gsw_cpp/gsw-node)
