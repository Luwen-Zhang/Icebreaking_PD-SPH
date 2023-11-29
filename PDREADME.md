# 写在前面

请在了解和学习了SPHinXsys之后，再参考以下修改记录

# 影响函数及其影响半径

..\SPHinXsys-master-source\SPHINXsys\src\shared\kernels\kernel_wenland_c2.h
从class Kernel派生出class Winfunc用以引入PD影响函数，kernel_size_保持为2.0
注意，增算影响函数的导数，用于计算固体表面法向量

..\SPHinXsys-master-source\SPHINXsys\src\shared\bodies\solid_body.h
从class RealBody派生出class PDBody，通过defineAdaptation<SPHAdaptation>(1.5075);实现PD独有的影响半径

..\SPHinXsys-master-source\SPHINXsys\src\shared\adaptations\adaptation.h
增设std::string body_name_;将其前6个字符“PDBody”作为构建class Winfunc的唯一标识符

..\SPHinXsys-master-source\SPHINXsys\src\shared\adaptations\adaptation.cpp
在class SPHAdaptation的构造函数中传参body_name_ = sph_body.getName();
并根据标识符“PDBody”，构建class Winfunc，
if (body_name_ == "PDBody") {
			kernel_ptr_ = makeUnique<Winfunc>(h_ref_);
		}

# 数据容器与数据结构

..\SPHinXsys-master-source\SPHINXsys\src\shared\particles\solid_particles.h
声明了class PDParticles : public ElasticSolidParticles这里添加了particleLive，bond，damage。专属PD的成员变量
声明了class NosbPDParticles : public PDParticles这里为NosbPD添加专属变量
声明了class NosbPDPlasticParticles : public NosbPDParticles这里为NosbPD弹塑性模型添加专属变量
..\SPHinXsys-master-source\SPHINXsys\src\shared\particles\solid_particles.cpp
定义了以上class的构造函数，初始化变量及其必要成员函数

..\SPHinXsys-master-source\SPHINXsys\src\shared\particle_neighborhood\neighborhood.h
增设了成员变量StdLargeVec<bool> bondLive_;以引入键的概念

..\SPHinXsys-master-source\SPHINXsys\src\shared\particle_neighborhood\neighborhood.cpp
针对class NeighborBuilder::createNeighbor和NeighborBuilder::initializeNeighbor的常量影响半径的成员函数
增加针对bondLive_的初始化程序

# 材料的本构关系

..\SPHinXsys-master-source\SPHINXsys\src\shared\materials\elastic_solid.h-cpp
声明了class HughesWingetSolid : public LinearElasticSolid实现小变形 大转动条件下的增量型本构
在cpp文件中定义了相应的构造函数和成员函数

..\SPHinXsys-master-source\SPHINXsys\src\shared\materials\inelastic_solid.h
声明了class PlasticSolidforPD : public HughesWingetSolid弹塑性本构关系的基类
通过设计成员变量，引入屈服应力、各项同性硬化模量及其内变量、随动硬化模量及其内变量
通过设计成员函数，引入屈服函数、各项同性硬化、随动硬化、弹塑性本构关系，这里均为虚函数，为接下来引入具体的NosbPD弹塑性本构做准备

声明了class J2PlasticityforPD : public PlasticSolidforPD基于应变加法分解和增量客观性的，考虑各项同性硬化和随动硬化的J2本构，采用关联的流动法则，能够应用于金属的模拟，注意仅适用于3D

声明了class DruckerPragerPlasticityforPD : public PlasticSolidforPD基于应变加法分解和增量客观性的，考虑各项同性硬化的DP本构，采用关联的流动法则，能够应用于土壤、岩石的模拟，在此主要应用于冰的模拟，注意仅适用于3D

以上均在相应cpp文件中定义了相应的构造函数和成员函数

# 计算方法与过程

## 数据容器

..\SPHinXsys-master-source\SPHINXsys\src\shared\particle_dynamics\solid_dynamics\elastic_dynamics.h
定义了typedef DataDelegateSimple<NosbPDParticles> NosbPDSolidDataSimple;
定义了typedef DataDelegateInner<NosbPDParticles> NosbPDSolidDataInner;
定义了typedef DataDelegateInner<NosbPDPlasticParticles> NosbPDPlasticSolidDataInner;

## 空间积分

..\SPHinXsys-master-source\SPHINXsys\src\shared\particle_dynamics\solid_dynamics\elastic_dynamics.h
声明了class PDTimeStepInitialization，并在相应cpp文件中重载了成员函数；施加重力加速度
声明了class NosbPDShapeMatrix，并在相应cpp文件中重载了成员函数；计算形张量K
声明了class NosbPDSecondStep，并在相应cpp文件中重载了成员函数；计算变形梯度F、欧拉速率梯度G及其中间变量，进而计算柯西应力、PK1、力态T0
声明了class NosbPDSecondStepPlastic，并在相应cpp文件中重载了成员函数；class NosbPDSecondStep的弹塑性版本，调用弹塑性本构关系完成力态T0的计算
声明了class NosbPDThirdStep，并在相应cpp文件中重载了成员函数；NosbPD控制方程右边第一项内力的空间积分，通过加权求和实现
声明了class LittleWoodHourGlassControl，并在相应cpp文件中重载了成员函数；背景键力形式的沙漏控制
声明了class PairNumericalDampingforPD，并在相应cpp文件中重载了成员函数；将TLSPH中的Kelvin-Voigt形式的人工粘弹性阻尼引入NosbPD

## 断键准则

..\SPHinXsys-master-source\SPHINXsys\src\shared\particle_dynamics\solid_dynamics\elastic_dynamics.h
声明了class NosbPDCheckBondLive，并在相应cpp文件中重载了成员函数；断键准则的基类
声明了class BondBreakByPrinStress，并在相应cpp文件中重载了成员函数；基于最大主应力的断键准则
声明了class BondBreakBySigma1andSigma3，并在相应cpp文件中重载了成员函数；基于最大主应力和最大剪应力的断键准则
声明了class BondBreakByPlasticStrain，并在相应cpp文件中重载了成员函数；基于最大塑性应变的断键准则

## 时间积分

..\SPHinXsys-master-source\SPHINXsys\src\shared\particle_dynamics\solid_dynamics\elastic_dynamics.h
声明了class NosbPDFirstStep，并在相应cpp文件中重载了成员函数；NosbPD的Position-Verlet时间积分，前1 / 2步进
声明了class NosbPDFourthStep，并在相应cpp文件中重载了成员函数；NosbPD的Position-Verlet时间积分，后1 / 2步进
声明了class ADRFirstStep，并在相应cpp文件中重载了成员函数；自适应动态松弛第一步
声明了class ADRSecondStep，并在相应cpp文件中重载了成员函数；自适应动态松弛第二步
声明了class NosbPDFourthStepWithADR，并在相应cpp文件中重载了成员函数；应用自适应动态松弛的NosbPD的Position-Verlet时间积分，后1 / 2步进

## 与流固耦合相关的算法

..\SPHinXsys-master-source\SPHINXsys\src\shared\particle_neighborhood\neighborhood.cpp
针对class NeighborBuilderContact 的构造函数为PDBody和WaterB相互作用编写了专属的Kernel选择策略
这里需要添加修改的原因：在SPHinXsys的多分辨率流固耦合策略中，根据论文中的思想，当两相采用不同分辨率时，只需仅按照流体域的影响半径搜索流固接触对，即可实现多分辨率流固耦合。由于在TLSPH和WCSPH的耦合计算中，TLSPH的光滑长度比是1.1，WCSPH的光滑长度比是1.3，因此原本SPHinXsys是通过比较光滑长度的大小来确定半径较大的一相是WCSPH流体，进而选择影响半径。但是，在NOSB-PD和WCSPH的耦合计算中，NOSB-PD的光滑长度比是1.5075，WCSPH的光滑长度比是1.3，这时即便NOSB-PD的分辨率较低，粒子初始间距较小，其影响半径未必小于WCSPH，这时仍然按照上述策略选择影响半径，有可能选为NOSB-PD的影响半径，与实际的思想不符。因此，此处通过class名，确凿地判定，哪一相是WCSPH，并将其光滑长度选为搜索FSI接触对的影响半径

..\SPHinXsys-master-source\SPHINXsys\src\shared\particle_dynamics\fluid_dynamics\fluid_dynamics_complex.h
声明了class BaseIntegration1stHalfWithWallforPD，并在相应cpp文件中重载了成员函数；源项：压力梯度，密度耗散；实现PD专属的流固耦合算法
声明了class BaseIntegration2ndHalfWithWallforPD，并在相应cpp文件中重载了成员函数；源项：压力耗散，速度散度；实现PD专属的流固耦合算法

..\SPHinXsys-master-source\SPHINXsys\src\shared\particle_dynamics\solid_dynamics\fluid_structure_interaction.h
声明了class BasePressureForceAccelerationFromFluidforPD，并在相应cpp文件中重载了成员函数；实现PD专属的流固耦合算法；添加了体积粘度阻尼；具有监测附近流场压力的功能，用于后处理

## 前处理

..\SPHinXsys-master-source\SPHINXsys\src\shared\particle_dynamics\base_local_dynamics.h
声明了class BaseSPHBodyRotation用于对生成之后的SPHBody粒子进行刚体转动，通过正交张量Q在某一坐标面内进行转动

## 后处理

..\SPHinXsys-master-source\SPHINXsys\src\shared\particles\solid_particles_variable.h
定义了typedef DataDelegateSimple<NosbPDParticles> PDSolidDataSimple;
定义了typedef DataDelegateSimple<NosbPDPlasticParticles> PDPlasticSolidDataSimple;

声明了class VonMisesStressforPD用于输出Mises等效应力
声明了class VonMisesPlasticStrainforPD用于输出Mises等效塑性应变，注意仅适用于3D

..\SPHinXsys-master-source\SPHINXsys\src\for_2D_build\particles\solid_particles_supplementary.cpp
..\SPHinXsys-master-source\SPHINXsys\src\for_3D_build\particles\solid_particles_supplementary.cpp
分别在2D和3D的cpp文件中重载了相应的成员函数




































