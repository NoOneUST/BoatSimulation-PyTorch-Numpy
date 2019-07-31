# BoatSimulation
### An original and useful code for boat simulation

-----------------------------------------------------

####        '`Discrete3DHull_Pytorch.py`' is the pytorch version of the main file running on GPU. <br>
    
You can use it to initialize a 3D boat with an arbitary function and to add ballasts and masts to your boat. 
    
After this, you will be able to calculate the center of mass and center of buoyancy of your boat, with which you can obtain its waterline and righting moment. 
    
In the end, the program will tell you whether your boat is stable and give you some instinct results.<br>

####        '`Discrete3DHull_Pytorch.py`' is the numpy version of the main file running on CPU. <br>

####        '`plotAlgebraic3DHull.py`' is a program file with which you can visualize any algebric boat you want.<br>

####        `Enjoy` this program! <br>


####    `Final Work`:
![image](https://github.com/NoOneUST/BoatSimulation/blob/master/images/8.jpg)


####    `Boat Test`:
https://drive.google.com/open?id=1xxnOqMBaFyNcgnLHjJZ1gzvoHRIJBtvC


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
![image](https://github.com/NoOneUST/BoatSimulation/blob/master/images/7.png)<br>

####    `Output 1:`
    This is Loop:  0
    
    COM:  237 237 81
    
    Weight Before Ballast: 0.00252125458791852
    Weight After Ballast:  0.9486632347106934
    Weight After Mast:  1.0147582292556763
    
    COM:  237 237 79 
    
    Processed: 0.02%     Loss: -1     waterOffsetLowerBound: -1312     waterOffsetUpperBound: 1968
    Processed: 0.04%     Loss: -1.0     waterOffsetLowerBound: -1312     waterOffsetUpperBound: 328
    Processed: 0.06%     Loss: 0.9241608381271362     waterOffsetLowerBound: -492     waterOffsetUpperBound: 328
    Processed: 0.08%     Loss: -0.6471373438835144     waterOffsetLowerBound: -492     waterOffsetUpperBound: -82
    Processed: 0.1%     Loss: 0.9207621216773987     waterOffsetLowerBound: -287     waterOffsetUpperBound: -82
    Processed: 0.12%     Loss: 0.19799742102622986     waterOffsetLowerBound: -184     waterOffsetUpperBound: -82
    Processed: 0.13999999999999999%     Loss: -0.25148385763168335     waterOffsetLowerBound: -184     waterOffsetUpperBound: -133
    Processed: 0.16%     Loss: -0.03376365080475807     waterOffsetLowerBound: -184     waterOffsetUpperBound: -158
    Processed: 0.18%     Loss: 0.08185473829507828     waterOffsetLowerBound: -171     waterOffsetUpperBound: -158
    Processed: 0.2%     Loss: 0.01958693563938141     waterOffsetLowerBound: -164     waterOffsetUpperBound: -158
    
    SubmergedVolume:  0.0010105937253683805
    Weight of Displaced Water:  1.0075619441922754
    Weight of Boat:  1.0147582292556763
    
    Done 
    
    Loss=  -0.7091647014021873 % 
    Offset =  -161
    Buoyancy Torque under  120  degrees:  -0.0449344621260616  N * M
    
    COM:  [237, 237, 577]
    COB:  [237, 266, 565]

    The height of Boat is:  0.1067
    The Length of Boat is:  0.198
    The Width of Boat is:  0.1782
    
    This is Loop:  1
    
    COM:  237 237 81
    
    Weight Before Ballast: 0.00252125458791852
    Weight After Ballast:  0.9486632347106934
    Weight After Mast:  1.0147582292556763
    
    COM:  237 237 79 
    
    Processed: 0.02%     Loss: -1     waterOffsetLowerBound: -1312     waterOffsetUpperBound: 1968
    Processed: 0.04%     Loss: -1.0     waterOffsetLowerBound: -1312     waterOffsetUpperBound: 328
    Processed: 0.06%     Loss: 0.9241608381271362     waterOffsetLowerBound: -492     waterOffsetUpperBound: 328
    Processed: 0.08%     Loss: -0.9119965434074402     waterOffsetLowerBound: -492     waterOffsetUpperBound: -82
    Processed: 0.1%     Loss: 0.9241608381271362     waterOffsetLowerBound: -287     waterOffsetUpperBound: -82
    Processed: 0.12%     Loss: 0.429568886756897     waterOffsetLowerBound: -184     waterOffsetUpperBound: -82
    Processed: 0.13999999999999999%     Loss: -0.40110132098197937     waterOffsetLowerBound: -184     waterOffsetUpperBound: -133
    Processed: 0.16%     Loss: -0.012636375613510609     waterOffsetLowerBound: -184     waterOffsetUpperBound: -158
    Processed: 0.18%     Loss: 0.20884087681770325     waterOffsetLowerBound: -171     waterOffsetUpperBound: -158
    Processed: 0.2%     Loss: 0.08868617564439774     waterOffsetLowerBound: -164     waterOffsetUpperBound: -158
    Processed: 0.22%     Loss: 0.037850271910429     waterOffsetLowerBound: -161     waterOffsetUpperBound: -158
    
    SubmergedVolume:  0.0010219018440693617
    Weight of Displaced Water:  1.0188361385371536
    Weight of Boat:  1.0147582292556763
    
    Done 
    
    Loss=  0.40186038240790367 % 
    Offset =  -159
    Buoyancy Torque under  140  degrees:  -0.08158023033216084  N * M
    
    COM:  [237, 237, 577]
    COB:  [237, 261, 560]
    
    RightMomentList:  [-0.0449344621260616, -0.08158023033216084]
