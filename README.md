# GSW-node

Nodejs interface to Gibbs-SeaWater (GSW) Oceanographic Toolbox in C++.

Link to the homepage: [Teos-10.org](https://www.teos-10.org/)

## Requirments
- Node: 12
- cmake
- gcc

Unfortunatly this library is only compatable with Node versions up to 12. It may build on later versions but tested working on node 12. 

Currently this has only been tested on linux Ubuntu 20. 

## TeosNode
TODO: add link to npm package manager

If you clone this repo you can install using:
```bash 
$ yarn
```
Then test with:
```
$ yarn test
```

### Use
The use of gsw_z_from_p is shown in the below example, this is an example from https://www.teos-10.org/pubs/gsw/html/gsw_z_from_p.html

**Note** if you see the definition of gsw_z_from_p in the source code you will notice that there are two optional arguments, in the bindings of this function these defaults do not transfer over, in the following typescript code you can see that these defaults need to be added .
```c++
double gsw_z_from_p(double p, double lat, double geo_strf_dyn_height=0.0,
								  double sea_surface_geopotential=0.0);
```

```ts
import { TeosBase } from 'gsw-node'

const teosBase = new TeosBase()
const p = 10
const lat = 4
const expectedResult = -0.099445834469453 * 1.0e+002
const result = teosBase.gsw_z_from_p(p, lat, 0, 0)

```

## TeosCpp
This is a copy of the source code from Teos-10, downloadable from https://www.teos-10.org/software.htm

This has been converted from use with codeblocks build chain to use cmake, which may be of use to others. You can build this with the following steps
### Building
```bash
$ cd TeosCpp && mkdir build
$ cd build
$ cmake ..
$ make
```
Run very very basic the example code using:
```bash
./test-teos
```