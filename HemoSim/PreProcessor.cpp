#include "PreProcessor.h"
#include "Domain.h"
#include "ModelParam.h"
#include "Element.h"
#include <ctype.h>

void PreProcessor::run(int argc, char *argv[]){
  char buf[BUFSIZE], *bf, *infile;   /* buf: character string with the name of the input file ".in";
                                     "BUFSIZ = 1024" (defined in all the standard library
                                     routines) is the number of characters into the array "buf".
                                     infile: pointer to the string of the name of the input file ".in" */
  int  i,j,nvar,n;              /* nvar: number of the variables in the problem (p,u,hd) 
                                Ndoms: number of the domains (vessels), initially set to 1 */
  FILE    *fp;                  /* "fp" is declared as a file pointer, which points to a structure containing information
                                about the file to be read or written. In this case is the input file ".in". */

  /* "infile" is defined as a string with the name of the input file. */
  infile = argv[argc-1];

  sprintf(buf,"%s.in",  strtok(infile,".")); /* This "stdio.h" function writes the output "%s.in" into
                                             the string "buf", terminated with '\0'. The return count
                                             does not include the '\0'. "strtok" is a function from
                                             the library "string.h" that searches "infile" for tokens
                                             delimited by characters from ".". Hence "buf" is a
                                             character array containing the name of the input file
                                             "input.in".*/
  if(!(fp = fopen(buf,"r"))){   /* If the input file "%s.in" is not found, print an error message and
                                terminate executing the program. If found, open it to start reading it. */
    fprintf(stderr,"Error in command line: input file %s doesn't exist. \n",buf);
    exit(-1);  /* Standard library function that terminates program execution. It
               also calls "fclose" for each open output file in order to flush
               out any buffered output. */
  }

  /* Read the parameter list in the input file ".in". */
  ReadParams(fp);

  nvar = 3;  /* Number of variables of the problem (area, velocity and hd) -> "nvar" = 3.
             It is used by the function "ReadBC". */

  /* Find "Mesh" in the input file. */
  findSection("Mesh",buf,fp);

  /* In the input file ".in", if after the word "Mesh" it is written "Ndomains" or "ndomains",
  read the number of domains of the problem. */
  if(strstr(buf,"ndomains")||strstr(buf,"Ndomains")){  /* "strstr" is a "string.h" function that returns a
                                                       pointer to first occurrence of string "ndomains"
                                                       in "buf". */
    if((bf = (strchr(buf,'='))) == (char *) NULL){  /*If after the "if" statement the string "=" is not
                                                    found, terminate with the following error message. */
      fprintf(stderr,"Error in mesh definition: number of domains not specified. \n");
      exit(1);
    }
    sscanf(bf,"=%d",&ModelParam::Ndoms);
  }

  /* "omega" is initialized as a structure type "domain" pointer, which points to space for an array of
  "Ndoms" objects, each of size "sizeof(Domain)". */
  ModelParam::omega = (Domain *)calloc(ModelParam::Ndoms,sizeof(Domain));

  /* Set up the time integration order, "Je", as the maximum between 1 and the numerical value of
  "INTTYPE" read from the parameter list in the input file. */
  //g_Je = max(1,g_inttype);

  /* For each element of each domain i=0,...,Ndoms, read the parameters that defined each domain and store
  them in "omega[i].U" by means of "Setup".  The information read is:
  - elemenent identification "id"
  - polynomial order "L"
  - quadrature order "q"
  - space locations of both sides of the domain "x"
  It is worked out the
  - Jacobian of the elemenetal mapping onto the standard element "jac"
  - inverse of the jacobian "rx".
  The following vectors are declared
  - global numbering of the domain "Z"
  - solution in the quadrature points of the physical space "h"
  - solution in the projected space "hj"
  - material properties "Beta" in each quadrature point
  - initial area "Ao" in each quadrature point
  - pointer to the next element of the domain "next". */
  for(i = 0; i < ModelParam::Ndoms; ++i)
    ModelParam::omega[i].U = Setup(fp, &(ModelParam::omega[i].nel), ModelParam::omega[i], i);

  /* Read bc's information in the input file ".in". It is declarated and defined in "Setup.C". It is passed
  the number of domains "Ndoms", the number of variables ("nvar" = 2), the pointer "omega" to the
  structure "domain" containing information about all the elements of each domain, and the pointer "fp"
  to the input file ".in". */
  ReadBC(ModelParam::Ndoms,nvar,ModelParam::omega,fp);

  // if(g_ode_solver==2) // Bond Graph Model
  SortNode(ModelParam::omega);

  /* The structure element "Uf" is defined for each element of each domain as being
  identical to the structure element "U". */
  for(n = 0; n < ModelParam::Ndoms; ++n){
    ModelParam::omega[n].Uf = (Element **)calloc(ModelParam::oneD_Je,sizeof(Element *));
    for(i = 0 ; i < ModelParam::oneD_Je; ++i){
      ModelParam::omega[n].Uf[i] = CopyElmt(ModelParam::omega[n].U,ModelParam::omega[n].nel);
      ModelParam::omega[n].Uf[i]->h  = &(ModelParam::oneD_Ufh[i][n*ModelParam::oneD_q]);
      ModelParam::omega[n].Uf[i]->hj = &(ModelParam::oneD_Ufhj[i][n*ModelParam::oneD_L]);
    }
  }

  /* The structure element "Uf" is defined for each element of each domain as being
  identical to the structure element "U". */
  for(n = 0; n < ModelParam::Ndoms; ++n){
    ModelParam::omega[n].Uf = (Element **)calloc(ModelParam::oneD_Je,sizeof(Element *));
    for(i = 0 ; i < ModelParam::oneD_Je; ++i){
      ModelParam::omega[n].Uf[i] = CopyElmt(ModelParam::omega[n].U,ModelParam::omega[n].nel);
      ModelParam::omega[n].Uf[i]->h  = &(ModelParam::oneD_Ufh[i][n*ModelParam::oneD_q]);
      ModelParam::omega[n].Uf[i]->hj = &(ModelParam::oneD_Ufhj[i][n*ModelParam::oneD_L]);
    }
  }

  /* Define QGmax = max(QGmax, U->q) and LGmax = max(LGmax, U->L) taking into accound all the elements of
  all the domains. */
  library_init(ModelParam::omega,ModelParam::Ndoms);

  /* The structure elements "A" and "Af" are defined for each element of each domain as being
  identical to the structure element "U". */
  for(n = 0; n < ModelParam::Ndoms; ++n){
    ModelParam::omega[n].A = CopyElmt(ModelParam::omega[n].U,ModelParam::omega[n].nel);
    ModelParam::omega[n].A->h = &(ModelParam::oneD_Ah[n*ModelParam::oneD_q]);
    ModelParam::omega[n].Af = (Element **)calloc(ModelParam::oneD_Je,sizeof(Element *));
    for(i = 0 ; i < ModelParam::oneD_Je; ++i){
      ModelParam::omega[n].Af[i] = CopyElmt(ModelParam::omega[n].U,ModelParam::omega[n].nel);
      ModelParam::omega[n].Af[i]->h  = &(ModelParam::oneD_Afh[i][n*ModelParam::oneD_q]);
      ModelParam::omega[n].Af[i]->hj = &(ModelParam::oneD_Afhj[i][n*ModelParam::oneD_L]);
    }
  }

  for(n = 0; n < ModelParam::Ndoms; ++n) {
    ModelParam::omega[n].Hd = CopyElmt(ModelParam::omega[n].U, ModelParam::omega[n].nel);
    ModelParam::omega[n].Hd->h = &(ModelParam::oneD_Hdh[n*ModelParam::oneD_q]);
  }

  /* Initialize the oldAhis & oldUhis vector to TOL for each Domain */
  for(n = 0; n < ModelParam::Ndoms; ++n) {
    ModelParam::omega[n].oldAhis = dvector(0,ModelParam::omega[n].A[0].q-1);
    ModelParam::omega[n].oldUhis = dvector(0,ModelParam::omega[n].U[0].q-1);
    ModelParam::omega[n].oldHdhis = dvector(0,ModelParam::omega[n].Hd[0].q-1);
    dfill(ModelParam::omega[n].A[0].q, ModelParam::riemannTol, ModelParam::omega[n].oldAhis, 1);
    dfill(ModelParam::omega[n].U[0].q, ModelParam::riemannTol, ModelParam::omega[n].oldUhis, 1);
    dfill(ModelParam::omega[n].Hd[0].q, ModelParam::riemannTol, ModelParam::omega[n].oldHdhis, 1);
  }

  /* Read initial area and flow values in each quadrature point of each element of each domain. For each
  element, area is stored in the structure element "A" and flow is stored in the structure element "U". */
  ReadIC(ModelParam::Ndoms,nvar,ModelParam::omega,fp);

  /* If "History" can be found in the input file it means that history points are defined in some domains. */
  if(findSection("History",buf,fp)){
    int hispts, domid, ndomhis;
    fgets(buf,BUFSIZ,fp);
    sscanf(buf,"%d",&ndomhis); // Read the number of domains with history points "ndomhis".
    for(j = 0; j < ndomhis; ++j){ // For each domain with history points,
      fgets(buf,BUFSIZ,fp);
      sscanf(buf,"%d%d",&hispts,&domid); /* read the number of history points "hispts" and the number
                                         identifying the domain "domid". */
      domid = domid -1;
      ModelParam::omega[domid].hispts = hispts;  /* Store the number of history points "hispts" in the structure domain
                                     "omega". */
      ModelParam::omega[domid].hisptsx = dvector(0,ModelParam::omega[domid].hispts-1); /* Read the axial coordinates of the history
                                                               points of each domain and store it in the
                                                               structure domain "omega". */
      for(i=0; i < ModelParam::omega[domid].hispts; ++i)
        fscanf(fp,"%lf",ModelParam::omega[domid].hisptsx+i);
      fgets(buf,BUFSIZ,fp); // Finish off rest line
    }
  }

  /* If "Solution" can be found in the input file, read the solution and store it in "omega->solution". */
  if(findSection("Solution",buf,fp)){
    fgets(buf,BUFSIZ,fp);
    if((bf = strchr(buf,'=')) == (char *) NULL)
      error_msg("solution not specified correctly");
    while(isspace(*++bf));
    ModelParam::omega->solution = strdup(bf);
  }
}

/**
\brief
It reads the parameter list in the input file ".in"
*/
void PreProcessor::ReadParams(FILE *fp){
  int  n, NSTEPS; /* Number of parameters in the list */
  char buf[BUFSIZ], value[25], name[25], c;  /* buf: character string with a numerical figure
                                             indicating the number of parameters in
                                             the parameter list. BUFSIZ is the number
                                             of characters into the array buf.
                                             value: vector with the numerical figures of
                                             the values of the parameters in the list.
                                             name: vector with the name of the parameters in
                                             the list. */
  double dt;  /* dt: time step */

  /* Set the file position "fp" to the beginning of the file. */
  rewind (fp);

  fgets (buf, BUFSIZ, fp);  /* This function defined in the library "stdio.h" reads the "BUFSIZ-1"
                            characters into the array "buf", stopping if a newline is encountered;
                            the newline is included in the array which is terminated by '\0'. */
  if(sscanf(buf, "%d", &n) != 1) /*If the string "buf" is not a numerical figure, the program
                                 is terminated with a message error. */
  {fputs("Error in the parameter list: can't read # of parameters. ", stderr);exit(-1);}

  /* It is read each one of the "n" parameters contained in the list. */
  while (n--) {
    fgets (buf, BUFSIZ, fp);
    if(sscanf(buf, "%25s%25s", value, name) == 2){
      if(!strncmp(name,"DT",2))
        ModelParam::dt = atof(value);
      else if(!strncmp(name,"Rho",3))
        ModelParam::rho = atof(value);
      else if(!strncmp(name,"Alpha",5))
        ModelParam::oneD_alpha = atof(value);
      else if(!strncmp(name,"NSTEPS",6))
        ModelParam::nSteps = atoi(value);
      else if(!strncmp(name,"HISSTEP",7))
        ModelParam::hisSteps = atoi(value);
      else if(!strncmp(name,"RELTOL",6))
        ModelParam::relTol = atof(value);
      else if(!strncmp(name,"ABSTOL",6))
        ModelParam::absTol = atof(value);
      else if(!strncmp(name,"SHOWCFL",7))
        ModelParam::showCFL = atoi(value);
      else if(!strncmp(name,"STEPLAPSE",9))
        ModelParam::showStepLapse = atoi(value);
      else if(!strncmp(name,"TAPER",5))
        ModelParam::taperRate = atof(value);
      else if(!strncmp(name,"GAMMA",5)){
        if(atof(value)==1)
          ModelParam::gammaI = 1;
        else if(atof(value)==2)
          ModelParam::gammaII = 1;
      }
      else if(!strncmp(name,"scale_lamda",11)){
        ModelParam::scale_lamda = atof(value);
        ModelParam::nonDim = 1;
      }
      else if(!strncmp(name,"scale_u0",8)){
        ModelParam::scale_u0 = atof(value);
        ModelParam::nonDim = 1;
      }
      else if(!strncmp(name,"scale_r0",8)){
        ModelParam::scale_r0 = atof(value);
        ModelParam::nonDim = 1;
      }
      else if(!strncmp(name,"UpdVisc",7))
        ModelParam::updVisc = atoi(value);
      else if(!strncmp(name,"Riemann",7))
        ModelParam::riemann = atoi(value);
      else if(!strncmp(name,"Binout",6))
        ModelParam::outMode = atoi(value);
      else if(!strncmp(name,"ODESolver",9))
        ModelParam::solverType = atoi(value);
      else if(!strncmp(name,"BifurLoss",9))
        ModelParam::bifurloss = atoi(value);
    }
  }

  /* Different Riemann solver tolerance for dimensionalized and nondimensionalized settings. */
  if(ModelParam::nonDim)
    ModelParam::riemannTol = 1e-5;
  else
    ModelParam::riemannTol = 1e-15;

  return;
}

/**
\brief
It returns a char pointer which points to the first position after the string "name" in the file pointed
by "fp" . It also sets "fp" to point to this same position in the file.
*/
char *PreProcessor::findSection (char *name, char *buf, FILE *fp)
{
  char *p;
  while (p = fgets (buf, BUFSIZ, fp))
    if(strstr (p, name))
      break;
  return p;
}

/**
\brief
For each element of each domain i=0,...,Ndoms, read the parameters that define each domain and store
them in "omega[i].U".  The information read is:
- elemenent identification "id"
- polynomial order "L"
- quadrature order "q"
- space locations of both sides of the domain "x".

It is worked out the
- Jacobian of the elemenetal mapping onto the standard element "jac"
- inverse of the jacobian "rx".

The following vectors are declared
- global numbering of the domain "Z"
- solution in the quadrature points of the physical space "h"
- solution in the projected space "hj"
- material properties "Beta" in each quadrature point
- initial area "Ao" in each quadrature point
- pointer to the next element of the domain "next".
*/
Element *PreProcessor::Setup(FILE *fp, int *nel, Domain &omega_, int n){
  register int i,j,k;
  int cnt;
  char     buf[BUFSIZ], *st;  /* "buf" is a char variable where character strings are stored
                              to be used later, specially to be converted into numerical
                              values or to check the presence of a particular character string
                              in the input file ".in". */
  Element  *U;  /* Pointer to a structure "element". */
  double   beta = 0, Ao = 0, *x,*z,*w; /* "beta" and "Ao" are initially defined
                                       by the "d-parameters" read in the
                                       parameter list by the function
                                       "ReadParam".
                                       x: pointer to the axial coordinates.
                                       z: zeros of the quadrature points.
                                       w: weigths of the quadrature points.*/
  double PI = 4.0*atan(1.0); // Number pi

  /* The numerical figure representing the number of elements of the domain is read from the input file,
  which is pointed by "fp", and assigned to the char variable "buf". */
  fgets(buf,BUFSIZ,fp);
  /* The numerical figure of "buf" is converted to a numerical value and assigned to the int pointer
  "nel" (number of elements). */
  sscanf(buf,"%d",nel);

  U = (Element *)calloc(nel[0],sizeof(Element));  /* "U" is defined as a structure type element pointer,
                                                  which points to space for an array of "nel[0]"
                                                  objects, each of size "sizeof(Element)". */

  for(i = 0; i < nel[0]; ++i){ /* For each one of the elements of the domain, */
    if(i) U[i-1].next = U+i;  /* If there are two or more elements in the domain,
                              a pointer to the next element structure of the
                              domain is set up. */

    fgets(buf,BUFSIZ,fp);  /* The numerical figures representing "x_lower", "x_upper", "L" and "q" are
                           read from the input file, which is pointed by "fp", and assigned to the
                           char variable "buf". */
    if(strstr(buf,"=")==NULL)  /* If no "=" sign is found in the line (in the previously defined
                               variable "buf"), the numerical figures of "x_lower", "x_upper",
                               "L" and "q" are read from the char variable "buf" and assigned to
                               the respective members of the structure element pointed by "U". */
                               sscanf(buf,"%lf%lf%d%d",U[i].x,U[i].x+1,&U[i].L,&U[i].q);
    else{
      fprintf(stderr,"Error in mesh definition: possibly need 'Area' or 'Eh' in the heading of the domain. \n");
      exit(1);
    }
    U[i].id = i;  /* The element identification "id" is set up. */

#ifndef NODES
    /* 
    In Jordi's thesis, the relation is Qmin=L+1
    In Karniadakis's book, the relation is Qmin=L+2/3, which is equal to Qmin=L+2
    */
    U[i].q = max(U[i].q,U[i].L+1);
#endif
    ModelParam::oneD_L = U[i].L;
    ModelParam::oneD_q = U[i].q;

    if(ModelParam::oneD_Uh == NULL){
      ModelParam::oneD_Uh   = (double*)malloc(ModelParam::Ndoms*ModelParam::oneD_q*sizeof(double));
      ModelParam::oneD_Uhj  = (double*)malloc(ModelParam::Ndoms*ModelParam::oneD_L*sizeof(double));
      ModelParam::oneD_Ah   = (double*)malloc(ModelParam::Ndoms*ModelParam::oneD_q*sizeof(double));
      ModelParam::oneD_Ahj  = (double*)malloc(ModelParam::Ndoms*ModelParam::oneD_L*sizeof(double));
      ModelParam::oneD_Hdh  = (double*)malloc(ModelParam::Ndoms*ModelParam::oneD_q*sizeof(double));
      ModelParam::oneD_Hdhj = (double*)malloc(ModelParam::Ndoms*ModelParam::oneD_L*sizeof(double));
      ModelParam::oneD_Ufh  = (double**)malloc(ModelParam::oneD_Je*sizeof(double*));
      ModelParam::oneD_Ufhj = (double**)malloc(ModelParam::oneD_Je*sizeof(double*));
      ModelParam::oneD_Afh  = (double**)malloc(ModelParam::oneD_Je*sizeof(double*));
      ModelParam::oneD_Afhj = (double**)malloc(ModelParam::oneD_Je*sizeof(double*));

      for(int m=0; m<ModelParam::oneD_Je; m++)
      {
        ModelParam::oneD_Ufh[m]  = (double*)malloc(ModelParam::Ndoms*ModelParam::oneD_q*sizeof(double));
        ModelParam::oneD_Ufhj[m] = (double*)malloc(ModelParam::Ndoms*ModelParam::oneD_L*sizeof(double));
        ModelParam::oneD_Afh[m]  = (double*)malloc(ModelParam::Ndoms*ModelParam::oneD_q*sizeof(double));
        ModelParam::oneD_Afhj[m] = (double*)malloc(ModelParam::Ndoms*ModelParam::oneD_L*sizeof(double));
      }
    }

    /* Work out the Jacobian of the elemenetal mapping onto the standard element "jac" and the
    inverse of the jacobian "rx". */
    U[i].rx  = 2.0/(U[i].x[1] - U[i].x[0]);  /* rx = 2.0/(x_upper - x_lower) */
    U[i].jac = (U[i].x[1] - U[i].x[0])/2.0;  /* jac = (x_upper - x_lower)/2.0 */

    /* It is allocated memory for the following members of the structure "element". One position
    for each quadrature point. */
    U[i].Z    = ivector(0,U[i].L-1);
    U[i].beta = dvector(0,U[i].q-1);
    U[i].Ao   = dvector(0,U[i].q-1);
    U[i].Hdo  = dvector(0,U[i].q-1);

    U[i].h    = &(ModelParam::oneD_Uh[n*ModelParam::oneD_q]);
    U[i].hj   = &(ModelParam::oneD_Uhj[n*ModelParam::oneD_L]);

    x = dvector(0,U[i].q-1);  /* Vector with the axial coordinates of the "q" quadrature
                              points where the Area or Eh need to be evaluated. */
    getzw(U[i].q,&z,&w,'a'); /* Function defined in "Nektar/Hlib/src/Basis.C" that makes a link
                             list for zeros and weights for the Jacobi polynomial in [0,1]
                             at the gauss points. Note: alpha > -1 , beta > -1. The char 'a'
                             is the direction of zeros/weights, which can also be 'b' or 'c'. */
    for(j=0; j<U[i].q; j++)
      x[j] = 0.5*((1-z[j])*U[i].x[0]+(1+z[j])*U[i].x[1]);  /* The location of the quadrature points
                                                           is calculated. */

    /* Init the Area */
    fgets(buf,BUFSIZ,fp);
    if((st = strstr(buf,"="))!=NULL){  /* The presence of "=" is checked. */
      // dfill(U[i].q,atof(++st), U[i].Ao, 1);
      double initAo = atof(++st);

      for(int k=0; k<U[i].q; ++k){
        U[i].Ao[k] = (1-ModelParam::taperRate*k)*initAo;
        // printf("init Ao[k]=%f\n", U[i].Ao[k]);
      }
    }
    else{  /* If "=" not present the program terminates with the following error message. */
      fprintf(stderr,"Error in mesh definition: could not find string for 'Area'. \n");
      exit(1);
    }

    /* Init the viscosity */
    omega_.visc = dvector(0,U[i].q-1);
    fgets(buf, BUFSIZ, fp);
    if((st = strstr(buf, "=")) != NULL) {
      dfill(U[i].q, atof(++st), omega_.visc, 1);
      // printf("%f\n", omega_.visc);
    }

    /* Init the Hd */
    fgets(buf, BUFSIZ, fp);
    if((st = strstr(buf, "=")) != NULL) {
      dfill(U[i].q, atof(++st), U[0].Hdo, 1);
      // printf("%f\n", omega_.Hd);
    }

    /* Init gamma*/
    omega_.gamma = dvector(0,U[i].q-1);
    if(ModelParam::gammaI || ModelParam::gammaII) {  /* If viscoelastic wall model is active, the viscous modulus in each quadrature point is read and allocated.
                            Added by panqing. 2012-4-10 */
      fgets(buf, BUFSIZ, fp);
      if((st = strstr(buf, "=")) != NULL) {
        dfill(U[i].q, atof(++st), omega_.gamma, 1);
        // printf("%f\n", omega_.Hd);
      }
    }

    /* Init Eh */
    fgets(buf,BUFSIZ,fp);
    if((st = strstr(buf,"="))!=NULL){ /* The presence of "=" is checked. */
      // vector_def("x",++st);
      // vector_set(U[i].q,x,U[i].beta);
      dfill(U[i].q, atof(++st), U[i].beta, 1);
      for(k = 0; k < U[i].q; ++k){ // For each one of the quadrature points of the elemental region
        U[i].beta[k] = U[i].beta[k]*pow(PI,0.5)/0.75/U[i].Ao[k]; // Beta = 4/3*E*h*sqrt(PI)/Ao for each quadrature point
        //fprintf(stderr,"Beta = %lf \n", U[i].beta[k]);
      }
    }
    else{ /* If "=" not present the program terminates with the following error message. */
      fprintf(stderr,"Error in mesh definition: could not find string for Eh. \n");
      exit(1);
    }

    /* Init Boundary Hd */
    fgets(buf, BUFSIZ, fp);
    if((st = strstr(buf, "=")) != NULL) {
      omega_.BHd = atof(++st);
    }

    /* Init Boundary SO2 */
    fgets(buf, BUFSIZ, fp);
    if((st = strstr(buf, "=")) != NULL) {
      omega_.BSO2 = atof(++st);
    }

    /* Init Boundary Jm */
    fgets(buf, BUFSIZ, fp);
    if((st = strstr(buf, "=")) != NULL) {
      omega_.BJm = atof(++st);
    }

    /* Init Boundary Jc */
    fgets(buf, BUFSIZ, fp);
    if((st = strstr(buf, "=")) != NULL) {
      omega_.BJc = atof(++st);
    }

    /* Init Wall Thickness */
    fgets(buf, BUFSIZ, fp);
    if((st = strstr(buf, "=")) != NULL) {
      omega_.WallTh = atof(++st);
    }

    /* Init Measured Velocity */
    fgets(buf, BUFSIZ, fp);
    if((st = strstr(buf, "=")) != NULL) {
      omega_.MesVel = atof(++st);
    }

    free(x); /* If the initial area is variable the space pointed to by x (vectore with the quadrature points coordinates) is
             deallocated. This function is defined in stdlib.h */
  }

  double Len = 0;
  for(i=0;i<nel[0];++i)
    Len += U[i].jac*2;
  if (ModelParam::solverType==ModelParam::RLC_EXP || ModelParam::solverType==ModelParam::RLC_IMP 
    || ModelParam::solverType==ModelParam::RC_EXP || ModelParam::solverType==ModelParam::RC_IMP){
    // Fill the 0D model vars    
    omega_.ZeroD_R = 8*omega_.visc[0]*Len/(U[0].Ao[0]*U[0].Ao[0]/PI);
    omega_.ZeroD_C = 2*U[0].Ao[0]/sqrt(U[0].Ao[0])*Len/U[0].beta[0];
    omega_.ZeroD_I = ModelParam::rho*Len/U[0].Ao[0];

    // Init nodes
    for(i=0; i<4; ++i){
      omega_.nodes[0][i] = 0;
      omega_.nodes[1][i] = 0;
    }
  }
  
  if (ModelParam::solverType==ModelParam::SS){
    omega_.J = (U[0].Ao[0]*U[0].Ao[0]/ModelParam::PI)/(8*omega_.visc[0]*Len);
  }

  // Check for Ao area continuity
  for(i = 1; i < nel[0]; ++i)
    if(fabs((U[i-1].Ao[U[i-1].q-1] - U[i].Ao[0])) > 1e-15)
      fprintf(stderr,"Warning in mesh definition: Ao has a jump larger than 1e-15 between"
      " element %d and %d.\n",i,i+1);

  /* Set up global numbering */
  for(i = 0, cnt = 0; i < nel[0]; ++i){
    for(j = 0; j < U[i].L; ++j,++cnt)
      U[i].Z[j] = cnt;
  }

  return U;
}

/**
\brief
It reads the boundary conditions definition in the input file .in.
*/
void PreProcessor::ReadBC(int Ndoms, int nvar, Domain *omega, FILE *fp){
  register int i,n,k;
  int      tot=0,init;
  const    int nel = omega->nel;
  char     buf[BUFSIZ],*bf;
  Element  *U = omega->U;

  findSection("Boundary",buf,fp);

  for(n = 0; n < Ndoms; ++n){ /* For all the domains of the problem, */
    /* Allocate space for the BC type, value and character string of each domain. */
    omega[n].bctype    = (char *)   calloc(nvar*2+1,sizeof(char));
    omega[n].bcval     = (double *) calloc(nvar*2,sizeof(double));
    omega[n].bcstring  = (char **)  calloc(nvar*2,sizeof(char *));

    for(i = 0 ; i < nvar; ++i){  /* For all the variables of the problem (currently 2), B.C. at the
                                 lhs of the domain "n" are read and set up in the structure "omega". */
      fgets(buf,BUFSIZ,fp);
      sscanf(buf,"%1s%lf",omega[n].bctype+i, omega[n].bcval +i); /* The BC type and value are read for each variable of the problem. */
      switch(omega[n].bctype[i]){
      case 'd': case 'a': case 'u': case 'q': case 'p': case 'm': case 'h': /* For the BC types "d", "a", "u", "q", "p", and "m" (input W, area, velocity,
                                                                            flow rate, pressure, and mix absorbing-reflecting W) defined as functions, read the character string with
                                                                            the function of time that defines the BC. */
        fgets(buf,BUFSIZ,fp);
        if((bf = (strchr(buf,'='))) == (char *) NULL){  /* If B.C. defined as functions and no "="
                                                        string found terminate with an error
                                                        message. */
          fprintf(stderr,"Error B.C. definition: illegal B.C. function definition for input area, velocity, flow rate or pressure.\n");
          exit(1);
        }
        while(isspace(*++bf)); /* Read the input area and velocity function. */
        omega[n].bcstring[i] =  strdup(bf);
        break;
      case 'B': case 'C': case 'J': case 'I': /* For the BC types "B", "C", "J" and "I" (domain connected to a splitting flow bifurcation,
                                              a merging flow bifurcation, to another vessel by a union or to another vessel by a "ischemic" union),
                                              read and store information about the connectivity of the domain in the network. */
        if(!(omega[n].bifur)) omega[n].bifur = imatrix(0,1,0,1);
        sscanf(buf,"%*1s%d%d",omega[n].bifur[0],omega[n].bifur[0]+1);
        omega[n].bifur[0][0]--;
        omega[n].bifur[0][1]--;
        break;
      }
    }

    for(i = 0 ; i < nvar; ++i){  /* For all the variables of the problem (2 actually), B.C. at the
                                 rhs of the domain "n" are read and set up in the structure "omega". */
      fgets(buf,BUFSIZ,fp);
      sscanf(buf,"%1s%lf",omega[n].bctype+nvar+i, omega[n].bcval+nvar+i); /* The BC type and value are read for each variable of the problem. */

      switch(omega[n].bctype[nvar+i]){
      case 'd': case 'a': case 'u': case 'q': case 'p': case 'm': case 'h': /* /* For the BC types "d", "a", "u", "q", "p", and "m" (input W, area, velocity,
                                                                            flow rate, pressure, and mix absorbing-reflecting W) defined as functions, read the character string with the 
                                                                            function of time that defines the BC. */

        fgets(buf,BUFSIZ,fp);
        if((bf = (strchr(buf,'='))) == (char *) NULL){  /* If B.C. defined as functions and no "="
                                                        string found terminate with an error
                                                        message. */
          fprintf(stderr,"Error B.C. definition: illegal B.C. function definition for input area, velocity, flow rate or pressure.\n");
          exit(1);
        }
        while(isspace(*++bf));
        omega[n].bcstring[nvar + i] =  strdup(bf);
        break;
      case 'B': case 'C': case 'J': case 'I':  /* For the BC types "B", "C", "J" and "I" (domain connected to a splitting flow bifurcation,
                                               a merging flow bifurcation, to another vessel by a union or to another vessel by a "ischemic" union),
                                               read and store information about the connectivity of the domain in the network. */
        if(!(omega[n].bifur)) omega[n].bifur = imatrix(0,1,0,1);
        sscanf(buf,"%*1s%d%d",omega[n].bifur[1],omega[n].bifur[1]+1);
        omega[n].bifur[1][0]--;
        omega[n].bifur[1][1]--;
        break;
      }
    }
  }
}

/**
\brief
It returns an element with the same properties as the input element U.
*/
Element *PreProcessor::CopyElmt(Element *U,int nel){
  register int i;
  Element *newE;

  newE = (Element *)calloc(nel,sizeof(Element));

  memcpy(newE,U,nel*sizeof(Element));

  for(i = 0; i < nel; ++i){
    if(i) newE[i-1].next = newE+i;
    // newE[i].h  = dvector(0,newE[i].q-1);
    // newE[i].hj = dvector(0,newE[i].L-1);
  }

  return newE;
}

/*!
  \brief Sort the node number for the 0D model
  Fill the omega.nodes[0/1][0]
    某段血管起点(0)和终点(1)的节点编号
    范围: from 1 to ModelParam::Nnodes
  Fill the omega.nodes[0/1][1,2,3]
    某段血管起点(0)和终点(1)所连接的三段血管的血管段编号
    取值的意义: -1: inlets, -2: outlets, -3: junction
    others: id of segment that linked to the node
*/
void PreProcessor::SortNode(Domain *omega){
  int i,j,k,m;
  int temp,node_cnt=0;

  // store all the id of the connected segments
  for(i=0;i<ModelParam::Ndoms;++i){
    for(m=0;m<6;m=m+3){
      // One of the connected segment is the current investigated segment
      omega[i].nodes[m/3][1] = i;
      // The other two depends on the connection matrix
      switch (omega[i].bctype[m]){
        // If it's a boundary segment, node number is set to -1
        case 'u': case 'q': case 'p': case 'h':
          omega[i].nodes[m/3][2] = -1;
          omega[i].nodes[m/3][3] = -1;
          break;
        case 'B': case 'C':
          // nodes[n][1,2,3] saves the id of connected segments
          omega[i].nodes[m/3][2] = omega[i].bifur[m/3][0];
          omega[i].nodes[m/3][3] = omega[i].bifur[m/3][1];
          break;
        case 'J':
          omega[i].nodes[m/3][2] = omega[i].bifur[m/3][0];
          omega[i].nodes[m/3][3] = -3;
          break;
        case 'R': case 'T':
          omega[i].nodes[m/3][2] = -2;
          omega[i].nodes[m/3][3] = -2;
          break;
      }

      // bubble sort
      for(j=0;j<3;++j){
        for(k=1;k<3-j;++k){
          if(omega[i].nodes[m/3][k]>omega[i].nodes[m/3][k+1]){
            temp=omega[i].nodes[m/3][k];
            omega[i].nodes[m/3][k]=omega[i].nodes[m/3][k+1];
            omega[i].nodes[m/3][k+1]=temp;
          }
        }
      }
      // printf("%d\t%d\t%d\n", omega[i].nodes[m/3][1], omega[i].nodes[m/3][2], omega[i].nodes[m/3][3]);
    }

    // Scan the omega to check if the nodes are new
    // start node
    for(j=0;j<ModelParam::Ndoms;++j){
      if(j==i)
        continue;
      else if(omega[i].nodes[0][1] == -1){
        omega[i].nodes[0][0] = node_cnt++;  // If it's the inlet, the node number is set to -1
        break;
      }
      else if(omega[j].nodes[0][1]==omega[i].nodes[0][1] && omega[j].nodes[0][2]==omega[i].nodes[0][2] && omega[j].nodes[0][3]==omega[i].nodes[0][3]){
        omega[i].nodes[0][0]=omega[j].nodes[0][0];
        break;
      }
      else if(omega[j].nodes[1][1]==omega[i].nodes[0][1] && omega[j].nodes[1][2]==omega[i].nodes[0][2] && omega[j].nodes[1][3]==omega[i].nodes[0][3]){
        omega[i].nodes[0][0]=omega[j].nodes[1][0];
        break;
      }
    }
    if(j==ModelParam::Ndoms)
      omega[i].nodes[0][0]=node_cnt++;

    // end node
    for(j=0;j<ModelParam::Ndoms;++j){
      if(j==i)
        continue;
      else if(omega[i].nodes[1][1] == -2){
        omega[i].nodes[1][0] = -2;  // If it's the outlet, the node number is set to -1
        break;
      }
      else if(omega[j].nodes[0][1]==omega[i].nodes[1][1] && omega[j].nodes[0][2]==omega[i].nodes[1][2] && omega[j].nodes[0][3]==omega[i].nodes[1][3]){
        omega[i].nodes[1][0]=omega[j].nodes[0][0];
        break;
      }
      else if(omega[j].nodes[1][1]==omega[i].nodes[1][1] && omega[j].nodes[1][2]==omega[i].nodes[1][2] && omega[j].nodes[1][3]==omega[i].nodes[1][3]){
        omega[i].nodes[1][0]=omega[j].nodes[1][0];
        break;
      }
    }
    if(j==ModelParam::Ndoms)
      omega[i].nodes[1][0]=node_cnt++;
  }
  ModelParam::Nnodes = node_cnt;  // Number of nodes in the network
}

/*!
\brief
It reads the initial area and flow values in each quadrature point of each element of each domain.
The initial conditions can be either constant values or functions of the axial direction. There is the
possibility to define the IC's using the value of "Ao" defined previously. For each element, the area
is stored in the structure element "A" and the flow is stored in the structure element "U".
*/
void PreProcessor::ReadIC(int Ndoms, int nfields, Domain *omega, FILE *fp){
  int    i,n,d;
  double *x,*z,*w;
  char   buf[BUFSIZ],*s;
  Element *U, *Uin[MAXFIELDS]; /* MAXFIELDS defined in "oneDcode.h" as 2: area & velocity. */
  findSection("Initial condition",buf,fp);
  x = dvector(0,ModelParam::oneD_MaxQ-1);

  fgets(buf,BUFSIZ,fp);
  if(strstr(buf,"Restart")||strstr(buf,"RESTART")){

  }
  else{
    rewind(fp);
    findSection("Initial condition",buf,fp);

    for(d = 0; d < Ndoms; ++d){ /* For all the domains, */
      Uin[0] = omega[d].A;
      Uin[1] = omega[d].U;
      Uin[2] = omega[d].Hd;

      for(n = 0; n < nfields; ++n){  /* In each domain, for all its variables, while "=" is found, the
                                     initial conditions are read and set up. */
        if(!(s = strchr(fgets(buf,BUFSIZ,fp),'=')))
        {fprintf(stderr,"Error in initial conditions: the following function is invalid. %s\n",s);
        exit(-1);}
        while (isspace(*++s));

        /* The initial conditions are set up for all the quadrature points along the domain. */
        for(U=Uin[n];U;U = U->next){
          getzw(U->q,&z,&w,'a');

          if(strstr(s,"Ao")) /* Use the Ao defined in the mesh generation as initial condition. The 
                             initial conditions can be either constant values or functions of the
                             axial direction. */
                             for(i = 0; i < U->q; ++i)
                               U->h[i] = U->Ao[i];
          else{
            for(i = 0; i < U->q; ++i)
              x[i] = 0.5*(1-z[i])*U->x[0] + 0.5*(1+z[i])*U->x[1];
            // vector_def("x",s);
            // vector_set(U->q,x,U->h);
            dfill(U->q, atof(s), U->h, 1);
          }	
        }
      }
    }
  }
  free(x);
}
