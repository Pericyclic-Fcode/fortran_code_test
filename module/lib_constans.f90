module lib_constans
        implicit none
        integer,private,parameter::sp=selected_real_kind(p=6,r=37)!获得单精度实型类别号
        integer,private,parameter::dp=selected_real_kind(p=13,r=200)!获得双精度实型类别号，想改精度在这里调
        !均采用国际单位制，有上下标的常数用下划线_后接上下标名称表示
        !数学常数 
        real(kind=8),parameter::pi=3.14159265358979
        real(kind=8),parameter::e_0=2.71828182845904!自然对数的底直接出现在公式中较少，加个下标0，区别于常用的元电荷
        
        !通用常数
        real(kind=8),parameter::h=6.62606957E-34!普朗克常数
        real(kind=8),parameter::h_bar=1.054571726E-34!=h/2pi
        real(kind=8),parameter::c=299792458
        real(kind=8),parameter::epsilon_0=8.854187817E-12!真空介电常数
        real(kind=8),parameter::mu_0=1.256637061E-6!真空磁导率
        
        !电磁常数
        real(kind=8),parameter::mu_B=9.27400968E-24!玻尔磁子
        real(kind=8),parameter::e=1.602176565E-19!元电荷
        
        !原子与核常数
        real(kind=8),parameter::m_p=1.672621777E-27!质子质量
        real(kind=8),parameter::m_n=1.674927351E-27!中子质量
        real(kind=8),parameter::m_e=9.10938291E-31!电子质量
        real(kind=8),parameter::alpha=7.29735257E-3!精细结构常数
        real(kind=8),parameter::R=10973731.57!里德伯常数，理想气体常数都会背而且更常用的是玻尔兹曼常数，所以不用单独的R表示
        
        !物理化学常数
        real(kind=8),parameter::k_B=1.3806488E-23!玻尔兹曼常数
        real(kind=8),parameter::N_A=6.02214129E23!阿伏加德罗常数
        real(kind=8),parameter::R_0=8.3144621!理想气体常数
        real(kind=8),parameter::F_0=96485.3365!法拉第常数
        
end module
