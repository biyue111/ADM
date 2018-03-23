The pseudocode of project 2

## Node 0

```
1. Init: prime_num = 0; finished_process_number = 0;
2. Waiting for the result of other nodes
3. while finished_process_number < PROCESS_NUMBER{ // Get result from a process
   prime_num += process_result;
   finished_process_number++;
   }
4. return prime_num;
```

## Ohter Nodes

```
1. Find all prime in [3,sqrt(N)]
2. Init: Calcul lower_value and high_value, 
         Create list all_number = [lower_value, high_value] and mark all even number as composite number
         
```
         
