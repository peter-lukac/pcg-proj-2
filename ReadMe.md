# Paralelní výpočty na GPU Projekt č.2: Implementace OpenACC kódu

step1: 3/3  
step2: 2.5/3  
step3: 1/3  
step4: 3/3  
nbody.txt: 3/3  

Celkovo 12.5/15  

Poznámky:

Step 2: Bylo by rychlejší použit redukci, 
Step 3: Špatně počítá těžiště pro 16384 částic, V každé redukci alokujete paměť, to pak vede na synchronizaci GPU