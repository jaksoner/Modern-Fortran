[toc]

# Modern Fortran

## 1、Fortran简介

### 1.1 实例：海啸模拟器

#### 1.1.1 浅水方程

$$
\frac{\partial \textbf {u}}{\partial t} + \textbf u \cdot \nabla \textbf u = -g \nabla h 
{\partial h\over \partial t} = -\nabla \cdot(\textbf u(H+h))
$$

浅水方程。第一个是动量守恒方程，第二个是质量守恒方程。$\bold u$ 是2-d速度向量，$\bold g$ 是中立加速度，$h$ 是水深，$H$ 是未扰动的水深，$t$是时间。

---

## 2、最小的可运行程序

我们将从模拟背景流引起的水在空间中的运动开始，而不改变其形状。

### 2.1 编译和运行第一个程序

```fortran
program hello
    print *, 'hello world!'
end program hello
```

编译非常简单，只需将源文件传递给编译器，还可以指定output（）可执行文件的名称：

` gfortran hello.f90 -o hello`

在这个过程中包含两步

- 1、编译：产生.o文件
- 2、链接：二进制文件.o不是可执行文件。

上述命令可以表示为以下两个命令

`gfortran  -c  hello.f90`  编译选项-c表示只编译，不链接。

`gfortran  hello.o  -o  hello`

### 2.2 模拟物体的运动

在这个阶段，我们将只模拟由于背景流而导致的对象（或流体）的移动。控制方程如下：
$$
{\part h \over \part t} + c  {\part h \over \part x} = 0
$$
离散：
$$
{h_i^{n+1}-h_i^n \over \Delta t} + c{h_i^n-h_{i-1}^{n}\over\Delta x} =0\\
\\

h_i^{n+1} = h_i^n -c{\Delta t\over \Delta x}(h_i^n-h_{i-1}^{n})\\
$$
令$dh(i) = h(i)-h(i-1)$

### 2.3 实现该程序



```fortran
program tsunami 

    implicit none

    integer :: i, n

    integer,parameter :: grid_size = 0
    integer,parameter :: num_time_steps = 100

    
    real,parameter :: dt = 1, dx = 1, c = 1  ! dt时间步长[s]    dx网格步长[m]   c背景流速度[m/s]

    real :: h(grid_size), dh(grid_size)

    integer,parameter :: icenter = 25    
    real,parameter :: decay = 0.02

    if ( grid_size < 1 ) stop "grid_size must be > 0"
    if ( dt <= 0 ) stop "time step dt must be > 0"
    if ( dx <= 0 ) stop "grid spacing dx must be > 0"
    if ( c <= 0 ) stop "background flow speed c must be > 0"

    !   初始化水深
    do concurrent (i = 1:grid_size)
        h(i) = exp(-decay * (i - icenter) ** 2)
    end do

    print *, 0, h

    time_loop: do n = 1, num_time_steps
        dh(1) = h(1) - h(grid_size)	！ 周期性边界条件

        do i = 2, grid_size
            dh(i) = h(i) - h(i-1)
        end do

        do i = 1, grid_size
            h(i) = h(i) - c * dt / dx *dh(i)
        end do

        print *, n, h

    end do time_loop

end program tsunami
```

​	

```matlab
clc,clear
h = xlsread('E:\C++_code\fortran\modern_fortran\chapter2\a.xlsx');
x = 1:100;
%%
figure(1)
for i = 1:101
    area( x, h(i,2:end),'facecolor','blue')
    axis([1 100 0 1.2])
    grid on
    pause(0.1)
end
%%
figure(2)
for j = 1:25:76
    area(x, h(j,2:end),'facecolor','blue')
    title('time step 0,25,50,75')
    axis([1 100 0 1.2])
    ylabel 'height'
    hold on
    grid on
 
end
```

---

## 3、编写可重构代码

本章介绍函数和子程序，引入另一个方程，计算水波的速度，重构代码

```fortran
program tsunami

  implicit none

  integer :: n
  integer, parameter :: grid_size = 100
  integer, parameter :: num_time_steps = 100

  real, parameter :: dt = 1, dx = 1, c = 1
  real :: h(grid_size)

  integer, parameter :: icenter = 25
  real, parameter :: decay = 0.02

  if (grid_size < 1) stop 'grid_size must be > 0'        
  if (dt <= 0) stop 'time step dt must be > 0'           
  if (dx <= 0) stop 'grid spacing dx must be > 0'        
  if (c <= 0) stop 'background flow speed c must be > 0' 

  call set_gaussian(h, icenter, decay) 
  print *, 0, h

  time_loop: do n = 1, num_time_steps
    h = h - c * diff(h) / dx * dt 
    print *, n, h
  end do time_loop

contains
  pure function diff(x) result(dx) 
    real, intent(in) :: x(:)
    real :: dx(size(x))
    integer :: im
    im = size(x)
    dx(1) = x(1) - x(im)
    dx(2:im) = x(2:im) - x(1:im-1)
  end function diff

  pure subroutine set_gaussian(x, icenter, decay) 
    real, intent(in out) :: x(:)
    integer, intent(in) :: icenter
    real, intent(in) :: decay
    integer :: i
    do concurrent(i = 1:size(x))
      x(i) = exp(-decay * (i - icenter)**2)
    end do
  end subroutine set_gaussian
  
end program tsunami
```

---

## 4、模块介绍

从内置模块中导入指定参数，eg:

```fortran
use iso_fortran_env, only: int32, real32
```

差分模块：

```fortran
module mod_diff

  use iso_fortran_env, only: int32, real32
  implicit none

  private
  public :: diff_upwind, diff_centered

contains

  pure function diff_centered(x) result(dx)
    ! Returns a 2nd order centered difference of a 1-d array,
    ! with periodic boundary condition.
    real(real32), intent(in) :: x(:)
    real(real32) :: dx(size(x))
    integer(int32) :: im
    im = size(x)
    dx(1) = x(2) - x(im)
    dx(im) = x(1) - x(im-1)
    dx(2:im-1) = x(3:im) - x(1:im-2)
    dx = 0.5 * dx
  end function diff_centered


  pure function diff_upwind(x) result(dx)
    ! Returns a 1st-order upstream finite difference of a 1-d array,
    ! with periodic boundary condition.
    real(real32), intent(in) :: x(:)
    real(real32) :: dx(size(x))
    integer(int32) :: im
    im = size(x)
    dx(1) = x(1) - x(im)
    dx(2:im) = x(2:im) - x(1:im-1)
  end function diff_upwind

end module mod_diff
```

初值模块：

```fortran
module mod_initial
    use iso_fortran_env
    implicit none
    
contains
    pure subroutine set_gaussian(x, icenter, decay) 
        real, intent(in out) :: x(:)
        integer, intent(in) :: icenter
        real, intent(in) :: decay
        integer :: i
        do concurrent(i = 1:size(x))
            x(i) = exp(-decay * (i - icenter)**2)
        end do
    end subroutine set_gaussian

end module mod_initial

```

主程序：

```fortran
program tsunami

  ! This version solves the non-linear 1-d shallow water equations:
  !
  !     du/dt + u du/dx = - g dh/dx
  ! 
  !     du/dt + u du/dx = - g dh/dx
  !
  ! The initial conditions and the finite difference calculation 
  ! are abstracted as external procedures.

  use iso_fortran_env, only: int32, real32
  use mod_diff, only: diff => diff_centered
  use mod_initial, only: set_gaussian

  implicit none

  integer(int32) :: n

  integer(int32), parameter :: grid_size = 100 ! grid size in x
  integer(int32), parameter :: num_time_steps = 5000 ! number of time steps

  real(real32), parameter :: dt = 0.02 ! time step [s]
  real(real32), parameter :: dx = 1 ! grid spacing [m]
  real(real32), parameter :: g = 9.8 ! gravitational acceleration [m/s^2]
  real(real32), parameter :: hmean = 10 ! mean water depth [m]

  real(real32) :: h(grid_size), u(grid_size)

  integer(int32), parameter :: icenter = 25
  real(real32), parameter :: decay = 0.02

  character(*), parameter :: fmt = '(i0,*(1x,es15.8e2))'

  ! check input parameter values
  if (grid_size <= 0) stop 'grid_size must be > 0'
  if (dt <= 0) stop 'time step dt must be > 0'
  if (dx <= 0) stop 'grid spacing dx must be > 0'

  ! initialize water height to a Gaussian blob
  call set_gaussian(h, icenter, decay)

  ! initialize water velocity to zero
  u = 0

  ! write initial state to screen
  print fmt, 0, h

  time_loop: do n = 1, num_time_steps

    ! compute u at next time step
    u = u - (u * diff(u) + g * diff(h)) / dx * dt

    ! compute h at next time step
    h = h - diff(u * (hmean + h)) / dx * dt

    ! write current state to screen
    print fmt, n, h

  end do time_loop

end program tsunami
```

---

## 5、用数组分析时间序列数据

### 5.1 用Fortran数组分析股票价格

#### 5.1.1 此次练习的目的

1. **找到表现最好和最坏的股票**
2. **辨别出风险股票**
3. **什么时候去买进或卖出**

[股票数据下载](https://github.com/modern-fortran/stock-prices)

### 5.2 找到表现最好和最坏的股票

计算从2000到2018年不同股票的收益来估计股票的表现。

**动态数组**基本操作：

```fortran
real, allocatable :: a(:)  ! 数组声明
allocate(a(10))  ! 给声明的数组分配内存
deallocate(a)  ! 释放内存
allocated(a)  ! 返回占用内存状态
```

主程序：

```fortran
program stock_gain

    use mod_arrays, only: reverse
    use mod_io, only: read_stock

    implicit none

    character(len=4), allocatable :: symbols(:)
    character(len=:), allocatable :: time(:)
    real, allocatable :: open(:), high(:), low(:), close(:), adjclose(:), volume(:)
    integer :: i, im, n
    real :: gain

    symbols = ['AAPL', 'AMZN', 'CRAY', 'CSCO', 'HPQ ',&
                'IBM ', 'INTC', 'MSFT', 'NVDA', 'ORCL']

    do n = 1, size(symbols)

        call read_stock('./data/' // trim(symbols(n)) //  '.csv', time,&
        open, high, low, close, adjclose, volume)

        adjclose = reverse(adjclose)
        gain = (adjclose(size(adjclose)) - adjclose(1))

        if (n == 1) then
            print *, time(size(time)) // ' through ' // time(1)
            print *, 'Symbol, Gain (USD), Relative gain (%)'
            print *, '-------------------------------------'
        end if

        print *, symbols(n), gain, nint(gain / adjclose(1) * 100)

    end do

end program stock_gain
```

模块：

- mod_io.f90

```fortran
module mod_io

  ! A helper module for parsing stock price data in csv format.

  use mod_alloc, only: alloc

  implicit none

  private
  public :: read_stock, write_stock

contains

  integer function num_records(filename)
    ! Return the number of records (lines) of a text file.
    character(len=*), intent(in) :: filename
    integer :: fileunit
    open(newunit=fileunit, file=filename)
    num_records = 0
    do
      read(unit=fileunit, fmt=*, end=1)
      num_records = num_records + 1
    end do
    1 continue
    close(unit=fileunit)
  end function num_records

  subroutine read_stock(filename, time, open, high, low, close, adjclose, volume)
    ! Read daily stock prices from a csv file.
    character(len=*), intent(in) :: filename
    character(len=:), allocatable, intent(in out) :: time(:)
    real, allocatable, intent(in out) :: open(:), high(:), low(:),&
                                         close(:), adjclose(:), volume(:)
    integer :: fileunit, n, nm
    nm = num_records(filename) - 1
    if (allocated(time)) deallocate(time)
    allocate(character(len=10) :: time(nm))
    call alloc(open, nm)
    call alloc(high, nm)
    call alloc(low, nm)
    call alloc(close, nm)
    call alloc(adjclose, nm)
    call alloc(volume, nm)
   
    open(newunit=fileunit, file=filename) ! open CSV文件
    read(fileunit, fmt=*, end=1) ! 跳过第一行
    do n = 1, nm
      read(fileunit, fmt=*, end=1) time(n), open(n),&
        high(n), low(n), close(n), adjclose(n), volume(n)  ! 逐行读取数据并将它们存储在数组中
    end do
    1 close(fileunit) ! 完成后关闭文件
  end subroutine read_stock

  subroutine write_stock(filename, time, price, mvavg, mvstd)
    ! Write derived stock data to file.
    character(len=*), intent(in) :: filename
    character(len=:), allocatable, intent(in) :: time(:)
    real, intent(in) :: price(:), mvavg(:), mvstd(:)
    integer :: fileunit, n
    open(newunit=fileunit, file=filename)
    do n = 1, size(time)
      write(fileunit, fmt=*) time(n), price(n), mvavg(n), mvstd(n)
    end do
    close(fileunit)
  end subroutine write_stock 

end module mod_io
```

- mod_alloc.f90

```fortran
module mod_alloc

  implicit none

  private
  public :: alloc, free

contains

  subroutine alloc(a, n)
    real, allocatable, intent(in out) :: a(:)
    integer, intent(in) :: n
    integer :: stat
    character(len=100) :: errmsg
    if (allocated(a)) call free(a)
    allocate(a(n), stat=stat, errmsg=errmsg)
    if (stat > 0) error stop errmsg
  end subroutine alloc

  subroutine free(a)
    real, allocatable, intent(in out) :: a(:)
    integer :: stat
    character(len=100) :: errmsg
    if (.not. allocated(a)) return
    deallocate(a, stat=stat, errmsg=errmsg)
    if (stat > 0) error stop errmsg
  end subroutine free

end module mod_alloc
```

- mod_arrays.f90

```fortran
module mod_arrays

  ! Utility functions that operate on arrays.

  implicit none

  private
  public :: argsort, average, crossneg, crosspos, intdate, moving_average,&
            moving_std, reverse, std

contains

  pure function argsort(x) result(a)
    ! Returns indices that sort x from low to high.
    real, intent(in):: x(:)
    integer :: a(size(x))
    integer :: i, i0, tmp1
    real :: tmp2
    real :: xwork(size(x))
    a = [(real(i), i = 1, size(x))]
    xwork = x
    do i = 1, size(x) - 1
      i0 = minloc(xwork(i:), 1) + i - 1
      if (i0 /= i) then
        tmp2 = xwork(i)
        xwork(i) = xwork(i0)
        xwork(i0) = tmp2
        tmp1 = a(i)
        a(i) = a(i0)
        a(i0) = tmp1
      end if
    end do
  end function argsort

  pure real function average(x)
    ! Returns a average of x.
    real, intent(in) :: x(:)
    average = sum(x) / size(x)
  end function average

  pure function crossneg(x, w) result(res)
    ! Returns indices where input array x crosses its
    ! moving average with window w from positive to negative.
    real, intent(in) :: x(:)
    integer, intent(in) :: w
    integer, allocatable :: res(:)
    real, allocatable :: xavg(:)
    logical, allocatable :: greater(:), smaller(:)
    integer :: i
    res = [(i, i = 2, size(x))]
    xavg = moving_average(x, w)
    greater = x > xavg
    smaller = x < xavg
    res = pack(res, smaller(2:) .and. greater(:size(x)-1))
  end function crossneg

  pure function crosspos(x, w) result(res)
    ! Returns indices where input array x crosses its
    ! moving average with window w from negative to positive.
    real, intent(in) :: x(:)
    integer, intent(in) :: w
    integer, allocatable :: res(:)
    real, allocatable :: xavg(:)
    logical, allocatable :: greater(:), smaller(:)
    integer :: i
    res = [(i, i = 2, size(x))]
    xavg = moving_average(x, w)
    greater = x > xavg
    smaller = x < xavg
    res = pack(res, greater(2:) .and. smaller(:size(x)-1))
  end function crosspos

  pure elemental integer function intdate(t)
    ! Converts a time stamp in format YYYY-mm-dd to integer.
    character(len=10), intent(in) :: t
    character(len=8) :: str
    str = t(1:4) // t(6:7) // t(9:10)
    read(str, *) intdate
  end function intdate

  pure function moving_average(x, w) result(res)
    ! Returns the moving average of x with one-sided window w.
    real, intent(in) :: x(:)
    integer, intent(in) :: w
    real :: res(size(x))
    integer :: i, i1
    do i = 1, size(x)
      i1 = max(i-w, 1)
      res(i) = average(x(i1:i))
    end do 
  end function moving_average

  pure function moving_std(x, w) result(res)
    ! Returns the moving standard deviation of x with one-sided window w.
    real, intent(in) :: x(:)
    integer, intent(in) :: w
    real :: res(size(x))
    integer :: i, i1
    do i = 1, size(x)
      i1 = max(i-w, 1)
      res(i) = std(x(i1:i))
    end do 
  end function moving_std

  pure function reverse(x)
    ! Reverses the order of elements of x.
    real, intent(in) :: x(:)
    real, allocatable :: reverse(:)
    reverse = x(size(x):1:-1)
  end function reverse

  pure real function std(x)
    ! Returns the standard deviation of x.
    real, intent(in) :: x(:)
    std = sqrt(average((x - average(x))**2))
  end function std

end module mod_arrays
```

结果如下：

![image-20220808195654460](C:\Users\Administrator\AppData\Roaming\Typora\typora-user-images\image-20220808195654460.png)

### 5.3 辨别出风险股

```fortran
program stock_volatility

  use mod_arrays, only: average, std, moving_average, moving_std, reverse
  use mod_io, only: read_stock, write_stock

  implicit none

  character(len=4), allocatable :: symbols(:)
  character(len=:), allocatable :: time(:)
  real, allocatable :: open(:), high(:), low(:), close(:), adjclose(:), volume(:)
  integer :: i, im, n

  symbols = ['AAPL', 'AMZN', 'CRAY', 'CSCO', 'HPQ ',&
             'IBM ', 'INTC', 'MSFT', 'NVDA', 'ORCL']

  do n = 1, size(symbols)

    call read_stock('data/' // trim(symbols(n)) //  '.csv', time,&
      open, high, low, close, adjclose, volume)

    im = size(time)
    adjclose = reverse(adjclose)

    if (n == 1) then
      print *, time(im) // ' through ' // time(1)
      print *, 'Symbol, Average (USD), Volatility (USD), Relative Volatility (%)'
      print *, '----------------------------------------------------------------'
    end if

    print *, symbols(n), average(adjclose), std(adjclose),&
      nint(std(adjclose) / average(adjclose) * 100)

    time = time(im:1:-1)

    call write_stock(trim(symbols(n)) // '_volatility.txt', time, adjclose,&
      moving_average(adjclose, 30), moving_std(adjclose, 30))

  end do

end program stock_volatility

```

### 5.4 合适时间段买入卖出

