# Online Parallel geometric attacks on neural key exchange protocol

Within this code, we have implemented Parallel geometric attack on key exchange protocol which has 2 layer neural networks. According to protocol, there are two side, let us say Network 1 and Network 2. Every network has unique architecture, it has fixed 2 layers, and it has just one output. But how many inputs they have and the range of weights depends on our will. During key exchange protocol, both two network assign their weights random in the specified borders. In every iteration, these two network used the same random inputs (these are -1 or +1) and produce their output. Both side known other’s output. So both side try to adjust their weights according to hebbian learning algorithm to produce the same output. If both networks produce the same output in for instance 50 successive iteration, we can called those networks as synchronized.
 
During these two networks synchronization  process, we assume that attackers somehow listen the network and they know the inputs and both networks output. We assume that the attackers also know the network architecture as well. Out attackers also tries to get synchronize both network. But we have handicap that the target is not fixed. It always change. So according to geometric attacks algorithm we updated weight in some special cases. For instance, in any iteration if both networks output are different, then attackers do not update their own weights. Attacker update their own weight if both networks output and attackers own  in any iteration are the same. In last case, if network1 and network2’s output are same but not attacker’s output, then we needs to select the perceptron which minimize the difference. So we need to update just selected perceptron weight not all weights. 
We implement our code according to MPI library. It means our attackers can run on different process if you want. In this way, the researchers claim that due to do diversity of random initial weight, we have an opportunity to increase the possibility to get synchronized network.

All perceptron in our network has to have the same number of input. It means general number of inputs should be divisible by number of perceptron. If not then we increase the number of input to get tose number divisible.


## Compile
```
$ mpicc -std=c99 -o test Source.c
```
## Run

```
$ ./test
```

## Results


Number of inputs : 20 </br>
Total numbers of perceptrons: 6 </br>
Note, due to the input per perceptron should be integer, adjusted number of input is : 24 </br>
Positive boundaries of the weights: 4 </br>
Number of attackers: 4 </br>

Initial state of networks </br>
First Network </br>
[4, -3, -2, 4] </br>
[-1, -3, 3, -1] </br>
[2, -2, -2, 3] </br>
[-4, -1, -4, -1] </br>
[-2, 2, -1, 3] </br>
[-3, 4, 2, 3] </br>

Second Network
[-4, 1, -4, 1]
[-3, 2, -3, 0]
[4, -3, 3, -3]
[4, -2, 1, -3]
[-2, -3, -1, -2]
[-4, 4, 3, 0]

Output 1: -1
Output 2: 1

Final State of Network 1 and 2

First Network
[-1, -2, 1, -1]
[0, -2, 0, 1]
[1, 1, -1, -2]
[0, 2, -1, 0]
[1, 1, 2, 1]
[-1, 2, 1, -1]

Second Network
[-1, -2, 1, -1]
[-1, -2, 0, 2]
[1, 0, 0, -2]
[0, 2, -2, 1]
[1, 0, 2, 0]
[-2, 2, 1, 0]

Output 1: 1
Output 2: -1

 Network1 and Network2 are NOT SYNCHRONISED !!!!
Total iter: 1000000
---  List of syncronised attackers --------
None.




Number of inputs : 40
Total numbers of perceptrons: 4
Note, due to the input per perceptron should be integer, adjusted number of inpu
t is : 40
Positive boundaries of the weights: 4
Number of attackers: 4
Initial state of networks
First Network
[3, 1, -3, 0, 4, 2, -4, 2, -3, -4]
[2, 2, -4, -3, -4, 2, 2, 4, 1, -4]
[0, 0, 2, 4, -2, 3, 4, -4, -3, -1]
[-2, -3, 4, -4, -4, -1, -4, 0, 2, 2]

Second Network
[2, 0, 3, 0, -3, -1, 3, -4, -1, -2]
[3, 3, 1, 1, 1, 2, -4, 4, 1, 4]
[1, -2, -2, 4, -3, -4, 4, 4, -1, -4]
[-4, -4, -4, -4, 3, 0, 3, -2, 0, 4]

Output 1: -1
Output 2: -1

Final State of Network 1 and 2

First Network
[0, 0, -1, 1, -1, -1, 3, 3, 0, -2]
[1, 0, 2, 1, -1, -1, -2, -3, 3, 2]
[-2, 1, -1, 0, -2, -3, 2, 1, -1, 1]
[1, 0, 3, 1, -1, 2, 2, 1, -2, 0]

Second Network
[0, 0, -1, 1, -1, -1, 3, 3, 0, -2]
[1, 0, 2, 1, -1, -1, -2, -3, 3, 2]
[-2, 1, -1, 0, -2, -3, 2, 1, -1, 1]
[1, 0, 3, 1, -1, 2, 2, 1, -2, 0]

Output 1: 1
Output 2: 1

 Network1 and Network2 are SYNCHRONISED !!
Synchronised iter: 50
Total iter: 113921
Time Taken: 449.000 ms
---  List of syncronised attackers --------
None

## References

[1] I. Kanter, W. Kinzel, E. Kanter, Secure exchange of information by synchronization of neural
networks, Europhys. Lett., 141, 2002.

[2] A. Ruttor, W. Kinzel, R. Naeh, I. Kanter, Genetic Attack on Neural Cryptography, arXiv:condmat/0512022
v2, June 2006




