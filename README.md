# BoatSimulation
### An original and useful code for boat simulation

-----------------------------------------------------

####        '`Discrete3DHull_Pytorch.py`' is the pytorch version of the main file running on GPU. <br>

Only this version is now still in maintance, make sure to use this as if you have GPU.
    
You can use it to initialize a 3D boat with an arbitary function and to add ballasts and masts to your boat. 
    
After this, you will be able to calculate the center of mass and center of buoyancy of your boat, with which you can obtain its waterline and righting moment. 
    
In the end, the program will tell you whether your boat is stable and give you some instinct results.<br>

####        '`Discrete3DHull_Pytorch.py`' is the numpy version of the main file running on CPU. <br>

####        '`plotAlgebraic3DHull.py`' is a program file with which you can visualize any algebric boat you want.<br>

####        `Enjoy` this program! <br>


####    `Final Work`:
![image](https://github.com/NoOneUST/BoatSimulation/blob/master/images/8.jpg)


####    `Boat Test`:
The boat test video is here:
https://drive.google.com/open?id=1xxnOqMBaFyNcgnLHjJZ1gzvoHRIJBtvC


####    `3D Print`:
All the 3D print source files you need, including solidworks file, STL generating tools written with matlab(This program is modified on other people's work), gcode file that could be printed directly on Ultimaker, are here:
https://drive.google.com/open?id=1Zk6Com191bygU-sJWlvQdeYMoig3mTV4

####    `Feedback`:
I know very clearly that this program is still far from maturation. If you have any questions or find any bugs, please feedback to me via email `lwangcg@connect.ust.hk`.<br><br>


####    `Notice`:

It's not permitted for students in HKUST's QEA class to copy any part of this code directly before they submit their code, or it will be treated as cheating.<br><br>


####    `Example Outputs`:

![image](https://github.com/NoOneUST/BoatSimulation/blob/master/images/1.png)
![image](https://github.com/NoOneUST/BoatSimulation/blob/master/images/2.png)
![image](https://github.com/NoOneUST/BoatSimulation/blob/master/images/3.png)
![image](https://github.com/NoOneUST/BoatSimulation/blob/master/images/4.png)
![image](https://github.com/NoOneUST/BoatSimulation/blob/master/images/5.png)
![image](https://github.com/NoOneUST/BoatSimulation/blob/master/images/6.png)
![image](https://github.com/NoOneUST/BoatSimulation/blob/master/images/7.png)
![image](https://github.com/NoOneUST/BoatSimulation/blob/master/images/12.png)
![image](https://github.com/NoOneUST/BoatSimulation/blob/master/images/10.png)
![image](https://github.com/NoOneUST/BoatSimulation/blob/master/images/11.png)
<br>

####    `Output:`
    This is Loop:  0

    COM:  237 237 81

    Weight Before Ballast: 0.00252125458791852
    Weight After Ballast:  1.2640440464019775
    Weight After Mast:  1.3301390409469604

    COM:  237 237 76 

    Processed: 0.02%     Loss: -1     waterOffsetLowerBound: -1312     waterOffsetUpperBound: 1968
    Processed: 0.04%     Loss: -1.0     waterOffsetLowerBound: -1312     waterOffsetUpperBound: 328
    Processed: 0.06%     Loss: 0.4679352939128876     waterOffsetLowerBound: -492     waterOffsetUpperBound: 328
    Processed: 0.08%     Loss: -0.7308024168014526     waterOffsetLowerBound: -492     waterOffsetUpperBound: -82
    Processed: 0.1%     Loss: 0.4653424322605133     waterOffsetLowerBound: -287     waterOffsetUpperBound: -82
    Processed: 0.12%     Loss: -0.08605213463306427     waterOffsetLowerBound: -287     waterOffsetUpperBound: -184
    Processed: 0.13999999999999999%     Loss: 0.2516588270664215     waterOffsetLowerBound: -235     waterOffsetUpperBound: -184
    Processed: 0.16%     Loss: 0.08304774761199951     waterOffsetLowerBound: -209     waterOffsetUpperBound: -184

    SubmergedVolume:  0.001328026526607573
    Weight of Displaced Water:  1.3240424470277503
    Weight of Boat:  1.3301390409469604

    Done 

    Loss=  -0.4583431873470545 % 
    Offset =  -196
    Buoyancy Torque under  120  degrees:  0.012794669516966006  N * M

    COM:  [237, 237, 580]
    COB:  [237, 256, 568]

    It can recover

    The height of Boat is:  0.1067
    The Length of Boat is:  0.198
    The Width of Boat is:  0.1782

    This is Loop:  1

    COM:  237 237 81

    Weight Before Ballast: 0.00252125458791852
    Weight After Ballast:  1.2640440464019775
    Weight After Mast:  1.3301390409469604

    COM:  237 237 76 

    Processed: 0.02%     Loss: -1     waterOffsetLowerBound: -1312     waterOffsetUpperBound: 1968
    Processed: 0.04%     Loss: -1.0     waterOffsetLowerBound: -1312     waterOffsetUpperBound: 328
    Processed: 0.06%     Loss: 0.4679352939128876     waterOffsetLowerBound: -492     waterOffsetUpperBound: 328
    Processed: 0.08%     Loss: -0.9328624606132507     waterOffsetLowerBound: -492     waterOffsetUpperBound: -82
    Processed: 0.1%     Loss: 0.4679352939128876     waterOffsetLowerBound: -287     waterOffsetUpperBound: -82
    Processed: 0.12%     Loss: 0.09061288833618164     waterOffsetLowerBound: -184     waterOffsetUpperBound: -82
    Processed: 0.13999999999999999%     Loss: -0.5431023836135864     waterOffsetLowerBound: -184     waterOffsetUpperBound: -133
    Processed: 0.16%     Loss: -0.24674388766288757     waterOffsetLowerBound: -184     waterOffsetUpperBound: -158
    Processed: 0.18%     Loss: -0.07777969539165497     waterOffsetLowerBound: -184     waterOffsetUpperBound: -171

    SubmergedVolume:  0.0013352405512705445
    Weight of Displaced Water:  1.331234829616733
    Weight of Boat:  1.3301390409469604

    Done 

    Loss=  0.0823802431114018 % 
    Offset =  -177
    Buoyancy Torque under  140  degrees:  -0.017293643294573724  N * M

    COM:  [237, 237, 580]
    COB:  [237, 252, 564]

    RightMomentList:  [0.012794669516966006, -0.017293643294573724]
