#define PI 3.1415926535897932384

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void StartUI(void);
void EndUI(void);
void Simulation(void);


void Initialization(void);
void Algorithm(void);
int CheckStop(void);
void MemorySaveResult(void);

int TimeUpdate(int count);
int Update(int count);
int RungeKutta(int count);