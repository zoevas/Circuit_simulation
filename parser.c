#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "parser.h"



/************************************************************************************************************************
 * This function adds a new element-node into the linked list. Returns NULL if the element allready exists or a pointer	*
 * to the newly added element-node otherwise.										*
 ************************************************************************************************************************/
struct element_t *addElement(char *name, short terminal_count, char *terminal_a, char *terminal_b, char *terminal_c, char* terminal_d, double value, double L, double W, struct transient_source_t *transient_source){
	struct element_t *new;
	
#ifdef CHECK_FOR_DUPLICATE_ELEMENTS
	// check if the element allready exists in the linked list
	for(new=list; new!=NULL; new=new->next) 
		if(strcmp(new->name, name)==0) return NULL;
#endif

	// allocate memory for the new element
	new = (struct element_t*)malloc(sizeof(struct element_t));
	MALLOC_RETURN_CHECK(new);

	new->name = (char*)malloc((strlen(name)+1)*sizeof(char));
	MALLOC_RETURN_CHECK(new->name);
		
	new->terminal_a = (char*)malloc((strlen(terminal_a)+1)*sizeof(char));
	MALLOC_RETURN_CHECK(new->terminal_a);

	new->terminal_b = (char*)malloc((strlen(terminal_b)+1)*sizeof(char));
	MALLOC_RETURN_CHECK(new->terminal_b);

	if(terminal_count>2){
		new->terminal_c = (char*)malloc((strlen(terminal_c)+1)*sizeof(char));
		MALLOC_RETURN_CHECK(new->terminal_c);
		}

	if(terminal_count>3){
		new->terminal_d = (char*)malloc((strlen(terminal_d)+1)*sizeof(char));
		MALLOC_RETURN_CHECK(new->terminal_d);
		}

	// copy the data to the new element
	strcpy(new->name, name);
	new->terminal_count = terminal_count;
	strcpy(new->terminal_a, terminal_a);
	strcpy(new->terminal_b, terminal_b);
	if(terminal_count>2) strcpy(new->terminal_c, terminal_c);
	if(terminal_count>3) strcpy(new->terminal_d, terminal_d);
	new->value = value;
	new->L = L;
	new->W = W;
	new->transient_source = transient_source;

	// link the new element to the linked list
	new->next = list;
	list = new;

	// increase the elements' count
	elements_count++;

	// return the pointer to the newly added element
	return new;
	}




/************************************************************************************************************************
 * This function prints the content of the linked list in a format similar to out input file's. Notice that it prints	*
 * the elements in an exactly reverce order than they were originaly read.						*
 ************************************************************************************************************************/
void printList(){
	struct element_t *i;

	for(i=list; i!=NULL; i=i->next){
		if(i->terminal_count==2) printf("2-terminal element: %s %s %s %lf \n", i->name, i->terminal_a, i->terminal_b, i->value);
		if(i->terminal_count==3) printf("3-terminal element: %s %s %s %s  \n", i->name, i->terminal_a, i->terminal_b, i->terminal_c);
		if(i->terminal_count==4) printf("4-terminal element: %s %s %s %s %s %lf %lf \n", i->name, i->terminal_a, i->terminal_b, i->terminal_c, i->terminal_d, i->L, i->W);
		}
	}




/************************************************************************************************************************
 * This function prints the transient sources and the respective values.						*
 ************************************************************************************************************************/
void printTransientSources(){
	int i;
	struct element_t *elem;

	for(elem=list; elem!=NULL; elem=elem->next){
		if(elem->transient_source!=NULL){
			switch(elem->transient_source->type){
			default:
				printf("DEBUG OUTPUT: error in transient source printing\n");
				break;
			case TRANSIENT_TYPE_EXP:
				printf("(%s) EXP: %2.2lf %2.2lf %2.2lf %2.2lf %2.2lf %2.2lf\n", elem->name, elem->transient_source->un.transient_exp.i1, elem->transient_source->un.transient_exp.i2, elem->transient_source->un.transient_exp.td1, elem->transient_source->un.transient_exp.tc1, elem->transient_source->un.transient_exp.td2, elem->transient_source->un.transient_exp.tc2);
				break;		
			case TRANSIENT_TYPE_SIN:
				printf("(%s) SIN: %2.2lf %2.2lf %2.2lf %2.2lf %2.2lf %2.2lf\n", elem->name, elem->transient_source->un.transient_sin.i1, elem->transient_source->un.transient_sin.ia, elem->transient_source->un.transient_sin.fr, elem->transient_source->un.transient_sin.td, elem->transient_source->un.transient_sin.df, elem->transient_source->un.transient_sin.ph);
				break;
			case TRANSIENT_TYPE_PULSE:
				printf("(%s) PULSE: %2.2lf %2.2lf %2.2lf %2.2lf %2.2lf %2.2lf %2.2lf\n", elem->name, elem->transient_source->un.transient_pulse.i1, elem->transient_source->un.transient_pulse.i2, elem->transient_source->un.transient_pulse.td, elem->transient_source->un.transient_pulse.tr, elem->transient_source->un.transient_pulse.tf, elem->transient_source->un.transient_pulse.pw, elem->transient_source->un.transient_pulse.per);
				break;
			case TRANSIENT_TYPE_PWL:
				printf("(%s) PWL: ", elem->name);
				for(i=0; i<elem->transient_source->un.transient_pwl.n; i++) printf("(%2.2lf %2.2lf) ", elem->transient_source->un.transient_pwl.t[i], elem->transient_source->un.transient_pwl.i[i]);
				printf("\n");
				break;
				}
			}
		}
	}




/************************************************************************************************************************
 * This function prints all the .DC and .TRAN statements.								*
 ************************************************************************************************************************/
void printDotStatements(){
	int i, j;

	printf("Total %d of .DC statements:\n", dotDCctr);
	for(i=0; i<dotDCctr; i++){
		printf("\t.DC %s %2.2lf %2.2lf %2.2lf\n\t.PLOT", dotDC[i].element_name, dotDC[i].start_value, dotDC[i].end_value, dotDC[i].increment);
		for(j=0; j<dotDC[i].plot_node_ctr; j++) printf(" V(%s)", dotDC[i].plot_node_name[j]);
		printf("\n");
		}

	printf("Total %d of .TRAN statements:\n", dotTRANflag);
	if(dotTRANflag){
		printf("\t.TRAN %2.2lf %2.2lf\n\t.PLOT", dotTRAN.time_step, dotTRAN.fin_time);
		for(j=0; j<dotTRAN.plot_node_ctr; j++) printf(" V(%s)", dotTRAN.plot_node_name[j]);
		printf("\n");
		}
	}




/************************************************************************************************************************
 * This function prints all the .OPTIONS statements									*
 ************************************************************************************************************************/
void printOptions(){
	printf("OPTIONS: ");
	if(optionSPD) printf("SPD ");
	if(optionITER) printf("ITER ");
	if(optionSPARSE) printf("SPARSE ");
	if(optionMETHOD == METHOD_TR) printf("METHOD=TR ");
	else if(optionMETHOD == METHOD_BE) printf("METHOD=BE ");
	printf("ITOL=%2.3lf\n", optionITOL);
	}




/************************************************************************************************************************
 * This functions opens the file whose filename is specified by char *fname. Reads the whole file line by line and 	*
 * adds a new node to the linked list for every new element it encounters. Returns PARSE_NOT_OK if an error occurs	*
 * or PARSE_OK otherwise.												*
 ************************************************************************************************************************/
int parse(char *fname){
	FILE *f;
	char line[64*MAX_STRING_SIZE], element_name[MAX_STRING_SIZE], terminal_a[MAX_STRING_SIZE], terminal_b[MAX_STRING_SIZE], terminal_c[MAX_STRING_SIZE], terminal_d[MAX_STRING_SIZE];
	char tmp[MAX_STRING_SIZE];
				
	int line_ctr=1, i=0, j=0, terminal_count, flagPlot;
	double element_value, element_L, element_W;
	struct transient_source_t *transient_source;

	// initializes the linked list
	list = NULL;

	// initializes the transient_source pointer
	transient_source = NULL;

	// initializes the voltage sources count and the inductors count
	voltage_sources_count = inductors_count = 0;

	// initializes the circuit's elements' count
	elements_count = 0;

	// initializes all the options and "dot statements"
	dotDCctr = 0;
	optionSPD = 0;
	optionITER = 0;
	optionSPARSE = 0;
	optionMETHOD = METHOD_TR;	//default if TR
	optionITOL = 1e-3;	//default value for ITOL
	dotTRANflag = 0;
	flagPlot = NO_PLOT;
	
	// open the input file specified as a console parameter
	f = fopen(fname, "r");
	if(f==NULL){
		perror(NULL);
		return PARSE_NOT_OK;
		}
	
	// get the next line from the input file
	while( fgets( line, sizeof(line), f) ){

		// ignore line if comment
		if(line[0] == '*'){
			line_ctr++;
			continue;	
			}

		// because of case insensitivity, always convert to capital
		if(line[0]>='a' && line[0]<='z') line[0] = line[0] - 'a' + 'A';
	
		// eat up the spare spaces or tab characters
		for(i=2; i<strlen(line); i++){
			if((line[i-1]==' ' || line[i-1]=='\t') && (line[i]==' ' || line[i]=='\t')){
				for(j=i; j<strlen(line); j++) 
					line[j] = line[j+1];
				i--;
				}
			if(line[i-1]=='\t') line[i-1] = ' ';
			}

		// ignore empty lines
		if(strlen(line)<=1){
			line_ctr++;
			continue;
			}

		// check for a .DC statement
		if(!strncmp(line, ".DC", 3) || !strncmp(line, ".dc", 3)){
			if(dotDCctr>=MAX_DC_STATEMENTS){
				printf("reached max number of .DC statments\n");
				return PARSE_NOT_OK;
				}
			if(sscanf(&line[3], " %s %lf %lf %lf", dotDC[dotDCctr].element_name, &(dotDC[dotDCctr].start_value), &(dotDC[dotDCctr].end_value), &(dotDC[dotDCctr].increment))==4){
				dotDC[dotDCctr].plot_node_ctr = 0;
				flagPlot = PLOT_OF_DC;
				}
			else{
				printf("error in .DC statement in line %d\n", line_ctr);
				return PARSE_NOT_OK;
				}
			}

		// check for a .TRAN statement
		else if(!strncmp(line, ".TRAN", 5) || !strncmp(line, ".tran", 5)){
			if(sscanf(&line[5], " %lf %lf", &(dotTRAN.time_step), &(dotTRAN.fin_time))==2){
				dotTRAN.plot_node_ctr = 0;
				flagPlot = PLOT_OF_TRAN;
				}
			else{
				printf("error in .TRAN statement in line %d\n", line_ctr);
				return PARSE_NOT_OK;
				}
			}

		// check for a .PLOT or .PRINT statement
		else if(!strncmp(line, ".PLOT", 5) || !strncmp(line, ".PRINT", 6) || !strncmp(line, ".plot", 5) || !strncmp(line, ".print", 6)){
			i=0;
			while(line[i++]!=' ');
			if(flagPlot == NO_PLOT){
				printf("warning in .PLOT statement without previous declaration of a .DC or a .TRAN statement in line %d\n", line_ctr);
				return PARSE_NOT_OK;
				}
			else if(flagPlot == PLOT_OF_DC){
				while(i<strlen(line) && sscanf(&line[i], "%*[V]%*[(]%[^)]%*[)]%*[ ,\t,\n]", dotDC[dotDCctr].plot_node_name[dotDC[dotDCctr].plot_node_ctr])){
					i+= strlen(dotDC[dotDCctr].plot_node_name[dotDC[dotDCctr].plot_node_ctr])+4;
					dotDC[dotDCctr].plot_node_ctr++;
					if(dotDC[dotDCctr].plot_node_ctr == MAX_PLOT_PER_DC){
						printf("reached max number of .PLOT statments per .DC\n");
						return PARSE_NOT_OK;
						}		
					}
				dotDCctr++;
				flagPlot = NO_PLOT;				
				}
			else if(flagPlot == PLOT_OF_TRAN){
				while(i<strlen(line) && sscanf(&line[i], "%*[V]%*[(]%[^)]%*[)]%*[ ,\t,\n]", dotTRAN.plot_node_name[dotTRAN.plot_node_ctr])){
					i+= strlen(dotTRAN.plot_node_name[dotTRAN.plot_node_ctr])+4;
					dotTRAN.plot_node_ctr++;
					if(dotTRAN.plot_node_ctr == MAX_PLOT_PER_TRAN){
						printf("reached max number of .PLOT statments per .TRAN\n");
						return PARSE_NOT_OK;
						}
					}
				dotTRANflag = 1;
				flagPlot = NO_PLOT;				
				}
			}

		// check for a .OPTIONS statement
		else if(!strncmp(line, ".OP", 3) || !strncmp(line, ".op", 3)){
			i=3;
			if(i<strlen(line) && !strncmp(&line[i], "TIONS", 5) || !strncmp(&line[i], "tions", 5)) i+=5;
			while(sscanf(&line[i], " %s", tmp)>0){
				if(!strcmp(tmp, "SPD") || !strcmp(tmp, "spd")){
					optionSPD = 1;
					}
				else if(!strcmp(tmp, "ITER") || !strcmp(tmp, "iter")){
					optionITER = 1;
					}
				else if(!strcmp(tmp, "SPARSE") || !strcmp(tmp, "sparse")){
					optionSPARSE = 1;
					}
				else if(!strncmp(tmp, "METHOD=", 7) || !strncmp(tmp, "method=", 7)){
					if(!strncmp(&tmp[7], "TR", 2) || !strncmp(&tmp[7], "tr", 2)) optionMETHOD = METHOD_TR;
					else if(!strncmp(&tmp[7], "BE", 2) || !strncmp(&tmp[7], "be", 2)) optionMETHOD = METHOD_BE;
					}
				sscanf(tmp, "ITOL=%lf", &optionITOL);
				sscanf(tmp, "itol=%lf", &optionITOL);
				i+= strlen(tmp);
				}
			flagPlot = NO_PLOT;
			}

		// ignore any unknown dot statemenets
		else if(line[0]=='.'){
			flagPlot = NO_PLOT;
			}			

		// check if we have a V,I,R,C,L or D element and parse it
		else if((line[0]=='V' || line[0]=='I' || line[0]=='R' || line[0]=='C' || line[0]=='L' || line[0]=='D') && sscanf(line, "%s %s %s %lf", element_name, terminal_a, terminal_b, &element_value) == 4){
			terminal_count = 2;
			element_L = element_W = 0.0;
			if(line[0]=='V') voltage_sources_count++;
			else if(line[0]=='L') inductors_count++;
			if(terminal_a[0] == 0 && terminal_b[0] == 0){
				printf("both element nodes are 0 in line %d\n", line_ctr);
				return PARSE_NOT_OK;
				}
			if(line[0]=='V' || line[0]=='I'){
				i = strlen(element_name) + strlen(terminal_a) + strlen(terminal_b) + 3;
				while(line[i++]!=' ');
				if(i<strlen(line) && !strncmp(&line[i], "EXP", 3)){
					transient_source = (struct transient_source_t*)malloc(sizeof(struct transient_source_t));
					MALLOC_RETURN_CHECK(transient_source);
					transient_source->type = TRANSIENT_TYPE_EXP;
					if(!sscanf(&line[i+3], " (%lf %lf %lf %lf %lf %lf)", &transient_source->un.transient_exp.i1, &transient_source->un.transient_exp.i2, &transient_source->un.transient_exp.td1, &transient_source->un.transient_exp.tc1, &transient_source->un.transient_exp.td2, &transient_source->un.transient_exp.tc2)){
						printf("error in transient definition line %d\n", line_ctr);
						return PARSE_NOT_OK;
						}
					}
				else if(i<strlen(line) && !strncmp(&line[i], "SIN", 3)){
					transient_source = (struct transient_source_t*)malloc(sizeof(struct transient_source_t));
					MALLOC_RETURN_CHECK(transient_source);
					transient_source->type = TRANSIENT_TYPE_SIN;
					if(!sscanf(&line[i+3], " (%lf %lf %lf %lf %lf %lf)", &transient_source->un.transient_sin.i1, &transient_source->un.transient_sin.ia, &transient_source->un.transient_sin.fr, &transient_source->un.transient_sin.td, &transient_source->un.transient_sin.df, &transient_source->un.transient_sin.ph)){
						printf("error in transient definition line %d\n", line_ctr);
						return PARSE_NOT_OK;
						}
					}
				else if(i<strlen(line) && !strncmp(&line[i], "PULSE", 5)){
					transient_source = (struct transient_source_t*)malloc(sizeof(struct transient_source_t));
					MALLOC_RETURN_CHECK(transient_source);
					transient_source->type = TRANSIENT_TYPE_PULSE;
					if(!sscanf(&line[i+5], " (%lf %lf %lf %lf %lf %lf %lf)", &transient_source->un.transient_pulse.i1, &transient_source->un.transient_pulse.i2, &transient_source->un.transient_pulse.td, &transient_source->un.transient_pulse.tr, &transient_source->un.transient_pulse.tf, &transient_source->un.transient_pulse.pw, &transient_source->un.transient_pulse.per)){
						printf("error in transient definition line %d\n", line_ctr);
						return PARSE_NOT_OK;
						}
					}
				else if(i<strlen(line) && !strncmp(&line[i], "PWL", 3)){
					transient_source = (struct transient_source_t*)malloc(sizeof(struct transient_source_t));
					MALLOC_RETURN_CHECK(transient_source);
					transient_source->type = TRANSIENT_TYPE_PWL;
					transient_source->un.transient_pwl.n = 0;
					i+=3;
					if(line[i]==' ') i++;
					while(i<strlen(line) && sscanf(&line[i], "(%lf %lf)%*[ ,\t,\n]", &transient_source->un.transient_pwl.t[transient_source->un.transient_pwl.n], &transient_source->un.transient_pwl.i[transient_source->un.transient_pwl.n])){
						while(line[i++]!=' ');
						while(line[i++]!=' ');
						transient_source->un.transient_pwl.n++;
						if(transient_source->un.transient_pwl.n == MAX_PWL_VALUES){
							printf("reached max number of PWL values\n");
							return PARSE_NOT_OK;
							}
						}
					}
				else transient_source = NULL;
				}
			flagPlot = NO_PLOT;
			}		
			
		// check if we have a Q element and parse it
		else if(line[0]=='Q' && sscanf(line, "%s %s %s %s", element_name, terminal_a, terminal_b, terminal_c) == 4){
			terminal_count = 3;
			element_value = element_L = element_W = 0.0;
			flagPlot = 0;
			}
			
		// check if we have a M element and parse it
		else if(line[0]=='M' && sscanf(line, "%s %s %s %s %s %lf %lf", element_name, terminal_a, terminal_b, terminal_c, terminal_d, &element_L, &element_W) == 7){
			terminal_count = 4;
			element_value = 0.0;
			flagPlot = NO_PLOT;
			}

		else{
			printf("error in element definition line %d\n", line_ctr);	
			return PARSE_NOT_OK;
			}

		// add the element to the list
		if(line[0]!='.' && addElement(element_name, terminal_count, terminal_a, terminal_b, terminal_c, terminal_d, element_value, element_L, element_W, transient_source)==NULL)
			printf("multiple declaration of element %s in line %d\n", element_name, line_ctr);

		// put transient_source pointer back to null
		transient_source = NULL;
		
		line_ctr++;
		}

	fclose(f);

	return PARSE_OK;
	}

