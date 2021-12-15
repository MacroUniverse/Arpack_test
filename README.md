* Arpack 主页 https://www.caam.rice.edu/software/ARPACK/
* Arpack++ 主页 http://www.ime.unicamp.br/~chico/arpack++/
* 然而这里的代码已经年久失修了

* 现在直接可以用的是 Arpack++2.3 https://github.com/m-reuter/arpackpp/releases/tag/2.3.0
* apt 安装： `sudo apt update` `sudo apt install libarpack++2-dev`, 会自动安装依赖的包
* 还要安装 gfortran 就可以了

* 进入 example 目录， 看见 Makefile 就直接 `make` 就行
* pdf 里面的第一个简单例子在 `examples/areig/nonsym/simple.cc` 里面， 运行 `make simple` 即可
