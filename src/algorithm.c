#include "algorithm.h"

/*========================= Variable Declaration =========================*/
double Tstart	=	0;
double Tfinal	=	2;
double Ts		=	0.0005;

double	WSensor			= 2*PI*100;
double	zetaSensor		= 0.7;
double	zetaActuator	= 0.7;


double	WActuatorTable[3] = {2*PI*16.4,	2*PI*20.5,	2*PI*24.6};
double	phicTable[3] = {0,	PI/8,	PI/4};

double WActuator, phic;

char	OutFileName[30] = {""};

double	time[4001]		= {0, };
double	phi[4001]		= {0, };
double	phihat[4001]	= {0, };
double	phihatdot[4001] = {0, };
double	p[4001]			= {0, };
double	phat[4001]		= {0, };
double	phatdot[4001]	= {0, };
double	delta[4001]		= {0, };
double	deltac[4001]	= {0, };
double	deltadot[4001]	= {0, };

// Runge=Kutta Variable
double k1_phi, k1_phihat, k1_phihatdot, k1_p, k1_phat, k1_phatdot, k1_delta, k1_deltac, k1_deltadot;
double k2_phi, k2_phihat, k2_phihatdot, k2_p, k2_phat, k2_phatdot, k2_delta, k2_deltac, k2_deltadot;
double k3_phi, k3_phihat, k3_phihatdot, k3_p, k3_phat, k3_phatdot, k3_delta, k3_deltac, k3_deltadot;
double k4_phi, k4_phihat, k4_phihatdot, k4_p, k4_phat, k4_phatdot, k4_delta, k4_deltac, k4_deltadot;

// Gain
double	Kp, Kd;

double	Lp, Lphi, Ldelta;

// Loop variable
int count, i, j;


void StartUI(void)
{
	printf("==================================================================\n");
	printf("|                  Numerical analysis Final Project              |\n");
	printf("==================================================================\n");

	WActuator	= WActuatorTable[2];
	phic		= phicTable		[0];

	printf("Frequency of Motor: %lf\n", WActuator);
	printf("Command of rocket roll: %lf\n", phic);
	printf("\nStart of Program. Press ENTER to start\n");
	getchar();

};


void EndUI(void)
{
	printf("\nWriting finished. Press ENTER to finish\n");
	getchar();
};



void Simulation(void)
{
	Initialization();

	do
	{
		// Time Update, Calculate K, Update state and propagate to next step
		Algorithm();
	}while( CheckStop() );	// Compare count and final step

	MemorySaveResult();
};



void Initialization(void)
{
	Kp				= 27.0;
	Kd				= 0.003;

	Lp				= 3.2;
	Lphi			= 1200;
	Ldelta			= 16000;

	count	= 0;

	phi[0] = phic + PI/180;
	deltac[0] = ( ( phic - phihat[0] ) * Kp - phat[0] ) * Kd;

	p[0] = 0.0;
	delta[0] = 0.0;
};

void Algorithm(void)
{
	// Time Update
	TimeUpdate(count);

	// Calculate each k
	RungeKutta(count);

	// Update, Propagation
	Update(count);
	
	count++;
};


int TimeUpdate(int count)
{
	// Require pre-determined Sampling period
	time[count] = Ts * count;

	return count;
};


int RungeKutta(int count)
{
	k1_phi			= Ts * p[count];
	k1_phihat		= Ts * phihatdot[count];
	k1_phihatdot	= Ts * ( -2 * zetaSensor * WSensor * phihatdot[count]		-		WSensor * WSensor * phihat[count]		+		WSensor * WSensor * phi[count]			);

	k1_p			= Ts * ( -Lp * p[count]										+		Lphi * sin( 4 * phi[count] )			+		Ldelta * delta[count]					);
	k1_phat			= Ts * phatdot[count];
	k1_phatdot		= Ts * ( -2 * zetaSensor * WSensor * phatdot[count]			-		WSensor * WSensor * phat[count]			+	 	WSensor * WSensor * p[count]			);

	k1_delta		= Ts * deltadot[count];
	k1_deltadot		= Ts * ( -2 * zetaActuator * WActuator * deltadot[count]	-		WActuator * WActuator * delta[count]	+		WActuator * WActuator * deltac[count]	) ;
	k1_deltac		= Ts * ( -phihatdot[count] * Kp * Kd						-		phatdot[count] * Kd																		); 


	k2_phi			= Ts * ( p[count]											+		k1_p/2 );
	k2_phihat		= Ts * ( phihatdot[count]									+		k1_phihatdot/2 );
	k2_phihatdot	= Ts * ( -2 * zetaSensor * WSensor * ( phihatdot[count]		+		k1_phihatdot/2 )						-		WSensor * WSensor * ( phihat[count]		+		 k1_phihat/2 )			+		WSensor * WSensor * ( phi[count]		+		k1_phi/2 ) ) ;

	k2_p			= Ts * ( -Lp * ( p[count]		+		 k1_p/2 )			+		Lphi * sin( 4 * ( phi[count]	+	 k1_phi/2 ) )			+		Ldelta * ( delta[count] + k1_delta/2 ) );
	k2_phat			= Ts * ( phatdot[count]		+		k1_phatdot/2 );
	k2_phatdot		= Ts * ( -2 * zetaSensor * WSensor * ( phatdot[count] + k1_phatdot/2 )	-	WSensor * WSensor * ( phat[count] + k1_phat/2 )		+		WSensor * WSensor * ( p[count] + k1_p/2 ) );

	k2_delta		= Ts * ( deltadot[count]		+	k1_deltadot/2 );
	k2_deltadot		= Ts * ( -2 * zetaActuator * WActuator * ( deltadot[count]	+	k1_deltadot/2 )			-		WActuator * WActuator * ( delta[count]	+	k1_delta/2 )		+		WActuator * WActuator * ( deltac[count] + k1_deltac/2 ) ) ;
	k2_deltac		= Ts * ( -( phihatdot[count]	+	k1_phihatdot/2 ) * Kp * Kd		-		(phatdot[count]		+	k1_phatdot/2 ) * Kd ); 


	k3_phi			= Ts * ( p[count]											+		k2_p/2 );
	k3_phihat		= Ts * ( phihatdot[count]									+		k2_phihatdot/2 );
	k3_phihatdot	= Ts * ( -2 * zetaSensor * WSensor * ( phihatdot[count]		+		k2_phihatdot/2 )						-		WSensor * WSensor * ( phihat[count]		+		 k2_phihat/2 )			+		WSensor * WSensor * ( phi[count]		+		k2_phi/2 ) ) ;

	k3_p			= Ts * ( -Lp * ( p[count]		+		 k2_p/2 )			+		Lphi * sin( 4 * ( phi[count]	+	 k2_phi/2 ) )			+		Ldelta * ( delta[count] + k2_delta/2 ) );
	k3_phat			= Ts * ( phatdot[count]		+		k2_phatdot/2 );
	k3_phatdot		= Ts * ( -2 * zetaSensor * WSensor * ( phatdot[count] + k2_phatdot/2 )	-	WSensor * WSensor * ( phat[count] + k2_phat/2 )		+		WSensor * WSensor * ( p[count] + k2_p/2 ) );

	k3_delta		= Ts * ( deltadot[count]		+	k2_deltadot/2 );
	k3_deltadot		= Ts * ( -2 * zetaActuator * WActuator * ( deltadot[count]	+	k2_deltadot/2 )			-		WActuator * WActuator * ( delta[count]	+	k2_delta/2 )		+		WActuator * WActuator * ( deltac[count] + k2_deltac/2 ) ) ;
	k3_deltac		= Ts * ( -( phihatdot[count]	+	k2_phihatdot/2 ) * Kp * Kd		-		(phatdot[count]		+	k2_phatdot/2 ) * Kd ); 


	k4_phi			= Ts * ( p[count]											+		k3_p );
	k4_phihat		= Ts * ( phihatdot[count]									+		k3_phihatdot );
	k4_phihatdot	= Ts * ( -2 * zetaSensor * WSensor * ( phihatdot[count]		+		k3_phihatdot )						-		WSensor * WSensor * ( phihat[count]		+		 k3_phihat )			+		WSensor * WSensor * ( phi[count]		+		k3_phi ) ) ;

	k4_p			= Ts * ( -Lp * ( p[count]		+		 k3_p )			+		Lphi * sin( 4 * ( phi[count]	+	 k3_phi ) )			+		Ldelta * ( delta[count] + k3_delta ) );
	k4_phat			= Ts * ( phatdot[count]		+		k3_phatdot );
	k4_phatdot		= Ts * ( -2 * zetaSensor * WSensor * ( phatdot[count] + k3_phatdot )	-	WSensor * WSensor * ( phat[count] + k3_phat )		+		WSensor * WSensor * ( p[count] + k3_p ) );

	k4_delta		= Ts * ( deltadot[count]		+	k3_deltadot );
	k4_deltadot		= Ts * ( -2 * zetaActuator * WActuator * ( deltadot[count]	+	k3_deltadot )			-		WActuator * WActuator * ( delta[count]	+	k3_delta )		+		WActuator * WActuator * ( deltac[count] + k3_deltac ) ) ;
	k4_deltac		= Ts * ( -( phihatdot[count]	+	k3_phihatdot ) * Kp * Kd		-		(phatdot[count]		+	k3_phatdot ) * Kd ); 

	return count;
};


int Update(int count)
{
	phi[count+1]		=	phi[count]			+	k1_phi/6		+	k2_phi/3		+	k3_phi/3		+	k4_phi/6;
	phihat[count+1]		=	phihat[count]		+	k1_phihat/6		+	k2_phihat/3		+	k3_phihat/3		+	k4_phihat/6; 
	phihatdot[count+1] 	= 	phihatdot[count]	+	k1_phihatdot/6	+	k2_phihatdot/3	+	k3_phihatdot/3	+	k4_phihatdot/6;

	p[count+1]			= 	p[count]			+	k1_p/6			+	k2_p/3			+	k3_p/3			+	k4_p/6;
	phat[count+1]		= 	phat[count]			+	k1_phat/6		+	k2_phat/3		+	k3_phat/3		+	k4_phat/6;
	phatdot[count+1] 	= 	phatdot[count]		+	k1_phatdot/6	+	k2_phatdot/3	+	k3_phatdot/3	+	k4_phatdot/6;

	delta[count+1] 		= 	delta[count]		+	k1_delta/6		+	k2_delta/3		+	k3_delta/3		+	k4_delta/6;
	deltadot[count+1] 	= 	deltadot[count]		+	k1_deltadot/6	+	k2_deltadot/3	+	k3_deltadot/3	+	k4_deltadot/6;
	deltac[count+1] 	= 	( ( phic			-	phihat[count+1] )	* Kp			-	phat[count+1] ) * Kd;	

	return count;
};



void MemorySaveResult(void)
{
	FILE* pFile;
	int PRINTER = 0;

	sprintf(OutFileName, "%10.6lf_%10.6lf", WActuator, phic) ;
	pFile = fopen( strcat(OutFileName, "_.txt"), "w+t" ) ;

	for( PRINTER = 0; PRINTER <= Tfinal/Ts; PRINTER++ )
		fprintf( pFile, "%6.4f\t\t%12.8f\t\t%12.8f\t\t%12.8f\t\t%12.8f\t\t%12.8f\n", time[PRINTER], phic*180/PI, phi[PRINTER]*180/PI, p[PRINTER]*180/PI, deltac[PRINTER]*180/PI, delta[PRINTER]*180/PI ) ;

	fcloseall();
};


int CheckStop(void)
{
	if(count < 4001) return 1;
	else             return 0;
};