# IIP speech preprocessing

C library for speech preprocessing.

## Target Platform
* Ubuntu 16.04 x64 (or higher)
* Windows 10 x64
* Mac high sierra 10.13.5
* Embedded board(raspberry pi, zynx board etc)

## Requirement
* cmake 3.5 or higher  
* (OPTION) OpenBLAS[Linux]  
* (OPTION) IntelMKL[Windows]  

## Notice
* Set Options of CMAKE properly for you Envrionment
* We don't provide backend library, we assume your Environment contains backend library
* If you want to use memory pool, you need to call 
```C
init(); //begining of your main()
finit(); //at the end of your main()
```

## Installation
[Installation Guide](https://github.com/gogyzzz/iip_sph_pp/wiki/Install_Guide)

## Example
[Code Example](https://github.com/gogyzzz/iip_sph_pp/wiki/Examples)

## Reference
* [Intel MKL Reference](https://software.intel.com/en-us/mkl-developer-reference-c-overview)
* [Netlib.org](http://www.netlib.org/)
* [Open audio library study document](https://github.com/kooBH/OpenAudioLibraryStudy)

## LICENSE
```
===========================================================
           Copyright (c) 2018, __IIPLAB__
                All rights reserved.

 This Source Code Form is subject to the terms of
 the Mozilla Public License, v. 2.0. 
 If a copy of the MPL was not distributed with this file,
 You can obtain one at http://mozilla.org/MPL/2.0/.
===========================================================
```

