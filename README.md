# QR Encoder

一年多前写的简易的二维码生成工具。刚学C语言时写的很丑的代码。可以帮助理解二维码的编码原理。

 ## 介绍

这是我在2016年初自己编写了基于C语言的QR code生成工具，目前可以将一串数字或ASCII码转换为QR code。

## 平台

- Windows
- Linux（不推荐，我花了一天才把allegro库给编译好。中途缺各种东西）

## 编译依赖

- allegro图形库

  当时见识短浅，不会写Windows桌面应用，更不会用Qt图形库，在网上找了这个偏僻的图形库，支持得不广泛，造成现在各种不兼容。

  Windows平台下allegro-5.0.10的下载地址：

  http://cdn.allegro.cc/file/library/allegro/5.0.10/allegro-5.0.10-mingw-4.7.0.zip

- TDM-GCC 4.7(这是Codeblocks 13.12内置的GCC)

  不保证其他版本的GCC能够编译。

## 在Codeblocks内配置allegro库

压缩文件解压后，将allegro库的lib文件夹里的所有*.a文件加入到Linker所链接的库里面(Compiler Settings->Linker Settings->Link Libraries)，以及把allegro库的include文件夹加入到compiler的搜索路径里(Compiler Settings->Search Directories->Compiler)。

## 运行时依赖

- libstdc++-6.dll(位于MinGW\bin\，bin目录在MinGW的installation guide中被推荐加入Path环境变量)
- allegro-5.0.10-md.dll(位于allegro\bin\) 



> 有问题请联系：songchaow@outlook.com