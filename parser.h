#ifndef PARSER_H
#define PARSER_H

#define CHECK_FOR_DUPLICATE_ELEMENTS

#define PARSE_OK	0
#define PARSE_NOT_OK	-1

#define MALLOC_RETURN_CHECK(X)\
	if(X==NULL){\
		perror(NULL);\
		exit(1);\
		}

#define MAX_STRING_SIZE		32
#define MAX_DC_STATEMENTS	16
#define MAX_PLOT_PER_DC		32
#define MAX_PLOT_PER_TRAN	32
#define MAX_PWL_VALUES		8

#define TRANSIENT_TYPE_EXP	1
#define TRANSIENT_TYPE_SIN	2
#define TRANSIENT_TYPE_PULSE	3
#define TRANSIENT_TYPE_PWL	4

#define NO_PLOT			0
#define PLOT_OF_DC		1
#define PLOT_OF_TRAN		2

#define METHOD_TR		0
#define METHOD_BE		1


short optionSPD;	// is TRUE if there is an SPD option used
short optionITER;	// is TRUE if there is an ITER option used
short optionSPARSE;	// is TRUE if there is an SPARSE option used
short optionMETHOD;	// is METHOD_TR or METHOD_BE

double optionITOL;



// the array of stacts that holds the info for all the parsed .DC statements
struct dotDC_t{
	char element_name[MAX_STRING_SIZE];
	double start_value;
	double end_value;
	double increment;
	char plot_node_name[MAX_PLOT_PER_DC][MAX_STRING_SIZE];	// an array to all the node names of the corresponding .PLOT statement
	int plot_node_ctr;	// the count of the node names in the corresponding .PLOT statement
	}dotDC[MAX_DC_STATEMENTS];

int dotDCctr;	// the count of .DC statements


// holds the info for the parsed .TRAN statement
struct dotTRAN_t{
	char element_name[MAX_STRING_SIZE];
	double time_step;
	double fin_time;
	char plot_node_name[MAX_PLOT_PER_TRAN][MAX_STRING_SIZE];	// an array to all the node names of the corresponding .PLOT statement
	int plot_node_ctr;	// the count of the node names in the corresponding .PLOT statement
	}dotTRAN;

short dotTRANflag;	// turns true if a .TRAN statement exists



// holds the transient spec of a transient source
struct transient_source_t{
	short type;
	union{
		struct transient_exp_t{
			double i1;
			double i2;
			double td1;
			double tc1;
			double td2;
			double tc2;
			}transient_exp;

		struct transient_sin_t{
			double i1;
			double ia;
			double fr;
			double td;
			double df;
			double ph;
			}transient_sin;

		struct transient_pulse_t{
			double i1;
			double i2;
			double td;
			double tr;
			double tf;
			double pw;
			double per;
			}transient_pulse;

		struct transient_pwl_t{
			int n;
			double t[MAX_PWL_VALUES];
			double i[MAX_PWL_VALUES];
			}transient_pwl;
		}un;
	};


int voltage_sources_count;	// the voltage sources count
int inductors_count;		// the inductors' count (used in MNA)
int elements_count;		// the count of all the circuit's elements


// this struct represents a circuit element stored in a single node of the linked list
struct element_t{
	char *name;		// starts with one of [V,I,R,C,L,D,M,Q]
	short terminal_count;	// the number of terminals of the element
	char *terminal_a;	// the positive terminal
	char *terminal_b;	// the negative terminal
	char *terminal_c;	// used by BJT and MOS
	char *terminal_d;	// used by MOS
	
	double value;
	double L;		// used in MOS
	double W;		// used in MOS

	int mna_terminal_a;	// the positive terminal used by the Modified Nodal Analysis (sequential node IDs)
	int mna_terminal_b;	// the negative terminal used by the Modified Nodal Analysis (sequential node IDs)

	struct transient_source_t *transient_source;

	struct element_t *next;
	};

struct element_t *list;	// the linked list's head

struct element_t *addElement(char *name, short terminal_count, char *terminal_a, char *terminal_b, char *terminal_c, char* terminal_d, double value, double L, double W,  struct transient_source_t *transient_source);

int parse(char *fname);	// the main parsing functions

void printList();
void printTransientSources();
void printDotStatements();
void printOptions();

#endif

