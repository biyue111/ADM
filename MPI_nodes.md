# 网络
在VMWare或VirtualBox打开两台或以上虚拟机（每台虚拟机需要有相同的用户名称，这里使用guest）
使用Host-Only配置网络（Internal理论上也可以，但是需要手动设置IP）

设置每个node的主机名（感觉可以不改）
  > sudo gedit /etc/hostname

更改Host文件（不是一定的，只是为了后面方便）

  > sudo gedit /etc/hosts
  
  添加
  > 192.168.1.133    node1
  
  > 192.168.1.135    node2

# 设置ssh免密登陆
默认使用ssh登陆其他主机（或自己）是需要密码的

## Node1(Server)
* 生成 SSH 私钥对，一路回车就好
  > ssh-keygen -t rsa 


* 生成authorized_keys并发送给其他node（这里举例node2）
  > ssh guest@node2 mkdir -p .ssh
  
  > cat .ssh/id_rsa.pub | ssh guest@node2 'cat >> .ssh/authorized_keys'
  
  （这里node2也可以直接是IP地址[就是如果前面不改hosts文件的话]）

## node2
* 确认连接到node1 （不要忘记exit）
  > ssh node1
  
# 配置机器数量以及进程数
在.c所在文件的位置添加hosts文件
输入
  > node1
  
  > node2
  
（这里node1或node2也可以直接是IP地址[就是如果前面不改hosts文件的话]）
  
或指定每个节点的进程数
  > node1 slots=4
  
  > node2 slots=2
  
最后确认每个node中.c文件均在相同目录名下
在node1，.c文件所在的文件夹中运行
  > mpiexec -hostfile hosts -np 8 ./cpi

（-np 8 在hosts文件中没有指定进程数是需要）

