/*
 * Nome: SchellingsModelMPI.c
 * Scopo: Simula il modello di segregazione proposto da Schelling
 * Autore: Giulio Triggiani
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "mpi.h"

#define DEMO 0      // Permette di mostrare la correttezza con una matrice fissata (0: non usare la demo, 1: usa la demo)

/*** Impostazione per la matrice di agenti ***/
#define ROWS 10                       // Numero di righe della matrice
#define COLUMNS 10                    // Numero di colonne della matrice
#define AGENT_X 'X'                    // Agente 'X'
#define AGENT_O 'O'                    // Agente 'O'
#define EMPTY ' '                      // Casella vuota
#define X_PERCENTAGE 30                // Percentuale di agenti 'X' nella matrice
#define O_PERCENTAGE 30                // Percentuale di agenti 'O' nella matrice
#define SAT_PERCENTAGE 33.333          // Percentuale di soddisfazione di un agente
#define MAX_STEP 100                   // Massimo numero di iterazioni
/*** Fine delle impostazioni per la matrice di genti ***/

/*** Altre impostazioni per la matrice ***/
#define MASTER 0                                          // Rank del processo master
#define SEED 15                                           // Seme per l'assegnazione delle celle libere
#define BLUE(string) "\033[1;34m" string "\x1b[0m"        // Colora di blu
#define RED(string) "\033[1;31m" string "\x1b[0m"         // Colora di rosso
/*** Fine delle impostazioni ***/

/*** Strutture per gestire la matrice ***/
typedef struct voidCell {
    int row_index;
    int column_index;
} voidCell;

typedef struct moveAgent {
    int destination_row;
    int destination_column;
    char agent;
} moveAgent;
/*** Fine delle strutture ***/

/*** Firme delle funzioni ***/
int generate_matrix(char *, int, int);                                                   // Funzione per generare ed inizializzare la matrice
int subdivide_matrix(int, int *, int *, int *);                                          // Funzione per suddividere la matrice tra i processi
void exchange_rows(int, int, int, char *, MPI_Comm);                                     // Funzione per scambiare le righe di ogni processo con i propri vicini
int *calculate_move(int, int, int, int, char *, int *);                                  // Funzione per calcolare gli agenti da spostare
int is_satisfied(int, int, int, int, int, int, char *);                                  // Funzione per controllare se un agente è soddisfatto (1: soddisfatto; 0: non soddisfatto)
voidCell *calculate_local_void_cells(int, char *, int, int *);                           // Funzione per calcolare le celle vuote locali ad un processo
voidCell *assign_void_cells(int, int, int, voidCell *, int *, MPI_Datatype, int);        // Funzione per unire tutte le celle vuote dei processi e restituire quelle di destinazione per il processo i-esimo
void move(int, int, int, char *, int *, voidCell *, int, int *, int *, MPI_Datatype);    // Funzione per spostare gli agenti
void calculate_total_satisfaction(int, int, char *);                                     // Funzione per calcolare la soddisfazione finale di tutti gli agenti della matrice

void define_voidCell_type(MPI_Datatype *);                                               // Funzione per definire il tipo voidCell
void define_moveAgent_type(MPI_Datatype *);                                              // Funzione per definire il tipo moveAgent
int calculate_source(int, int *, int *, int);                                            // Funzione per calcolare a quale processo appartiene una determinata riga della matrice
void synchronize(int, int, int *, int, moveAgent **, int, char *, MPI_Datatype);         // Funzione per sincronizzare gli spostamenti tra i processi
void print_matrix(int, int, char *);                                                     // Funzione per stampare la matrice
void err_finish(int *, int *, int *);                                                    // Funzione per terminare l'esecuzione in caso di errori

// DEMO
void test_init_matrix(char *matrix, int O_pct, int X_pct);
// DEMO
/*** Fine delle firme ***/

/*** Funzione main ***/
int main(int argc, char **argv) {
    // Definizione variabili
    int world_size, rank;                   // Numero di processori e rank del processore
    double start_time, end_time;            // Tempo di inizio e fine computazione
    char *matrix = NULL;                    // Matrice di char ('X', 'O', ' ')
    char *sub_matrix = NULL;                // Sottomatrice assegnata ad un processo
    int *displacements = NULL;              // Array che contiene i displacements per ogni processo (per la MPI_Scatterv)
    int *sendcounts = NULL;                 // Array che contiene il numero di elementi (#righe_assegnate * #colonne) di un processo
    int *rows_per_process = NULL;           // Array che contiene il numero di righe assegnate ad ogni processo
    int *want_move = NULL;                  // Array che indica quali agenti della sottomatrice vogliono muoversi
    int unsatisfied_agents = 0;             // Numero di agenti insoddisfatti per ogni processo (ad ogni iterazione)
    int number_of_local_void_cells = 0;     // Numero di celle vuote nella sottomatrice
    voidCell *local_void_cells = NULL;      // Array che contiene le celle vuote della sottomatrice
    int number_of_destination_cells = 0;    // Numero di celle vuote che sono state assegnate al processo
    voidCell *destinations = NULL;          // Array che contiene le celle vuote che sono state assegnate al processo dove poter spostare gli agenti

    // Inizializzazione MPI
    MPI_Status status;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);         // Rank del processo chiamante nel gruppo
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);   // Numero di processi nel gruppo

    MPI_Barrier(MPI_COMM_WORLD);                  // *****Tutti i processi sono inizializzati*****

    // Inizzilizzazione variabili
    start_time = MPI_Wtime();                                 // Ritorna il tempo passato dalla chiamata di un processo
    sendcounts = calloc(world_size, sizeof(int));
    displacements = calloc(world_size, sizeof(int));
    rows_per_process = calloc(world_size, sizeof(int));

    // Definizione MPI_Datatype
    MPI_Datatype VOID_CELL_TYPE;
    define_voidCell_type(&VOID_CELL_TYPE);
    MPI_Datatype MOVE_AGENT_TYPE;
    define_moveAgent_type(&MOVE_AGENT_TYPE);

    // Inizializzazione matrice
    if (world_size <= ROWS) {
        if (rank == MASTER) {
            matrix = malloc(ROWS * COLUMNS * sizeof(char *));
            if (DEMO)
                test_init_matrix(matrix, O_PERCENTAGE, X_PERCENTAGE);
            else if (!generate_matrix(matrix, O_PERCENTAGE, X_PERCENTAGE))
                err_finish(sendcounts, displacements, rows_per_process);

            // Mostra informazioni
            time_t mytime;
            struct tm *timeinfo;
            time(&mytime);
            timeinfo = localtime(&mytime);
            printf("Started at: %s", asctime(timeinfo));  
            printf("Number of processes: %d\n", world_size);
            printf("Number of iterations: %d\n", MAX_STEP);
            printf("\nMatrice iniziale:\n");
            print_matrix(ROWS, COLUMNS, matrix);
        }
    }

    // Calcolo della porzione di matrice da assegnare a ciascun processo
    if (!subdivide_matrix(world_size, displacements, sendcounts, rows_per_process))
        err_finish(sendcounts, displacements, rows_per_process);

    // Suddivisione delle righe tra i processi
    sub_matrix = malloc(rows_per_process[rank] * COLUMNS * sizeof(char *));
    MPI_Scatterv(matrix, sendcounts, displacements, MPI_CHAR, sub_matrix, rows_per_process[rank] * COLUMNS, MPI_CHAR, MASTER, MPI_COMM_WORLD);    // Funzione MPI che permette di dividere il carico sui processi (sendbuf, sendcounts, displacements, sendtype, recvbuf, recvcount, recvtype, root, comm)

    // Calcolo di quante righe 'originali' ha il processo e di quante ne ha 'totali'
    int total_rows = rows_per_process[rank];      // Righe con gia assegnate quelle in più
    int original_rows = total_rows - ((rank == 0 || rank == world_size - 1) ? 1 : 2);

    // Comincia l'esecuzione (verrà eseguita un massimo di MAX_STEP volte)
    for (int i = 0; i < MAX_STEP; i++) {
        // Scambia le righe tra i processi vicini e calcola gli agenti che si vogliono spostare
        exchange_rows(rank, world_size, original_rows, sub_matrix, MPI_COMM_WORLD);
        want_move = calculate_move(rank, world_size, original_rows, total_rows, sub_matrix, &unsatisfied_agents);

        // Calcolo delle celle vuote di ogni processo e assegnazione delle celle vuote a ciascun processo
        local_void_cells = calculate_local_void_cells(original_rows, sub_matrix, displacements[rank], &number_of_local_void_cells);
        destinations = assign_void_cells(rank, world_size, number_of_local_void_cells, local_void_cells, &number_of_destination_cells, VOID_CELL_TYPE, unsatisfied_agents);

        // Gli agenti insoddisfatti vengono spostati
        move(rank, world_size, original_rows, sub_matrix, want_move, destinations, number_of_destination_cells, displacements, sendcounts, MOVE_AGENT_TYPE);

        MPI_Barrier(MPI_COMM_WORLD);

        free(want_move);
        free(local_void_cells);
        free(destinations);
    }

    // Si recupera la matrice finale
    MPI_Gatherv(sub_matrix, sendcounts[rank], MPI_CHAR, matrix, sendcounts, displacements, MPI_CHAR, MASTER, MPI_COMM_WORLD);  // (sendbuff, sendcount, datatype, destbuff, destcount, displacements, datatype, root, comm)

    end_time = MPI_Wtime();
    MPI_Type_free(&VOID_CELL_TYPE);
    MPI_Type_free(&MOVE_AGENT_TYPE);
    MPI_Finalize();

    // Stampa matrice finale e calcolo della soddisfazione totale
    if (rank == MASTER) {
        printf("\nMatrice finale:\n");
        print_matrix(ROWS, COLUMNS, matrix);
        calculate_total_satisfaction(rank, world_size, matrix);
        printf("Time in ms = %f\n", end_time - start_time);
    }

    free(matrix);
    free(sub_matrix);
    free(sendcounts);
    free(displacements);
    free(rows_per_process);

    return 0;
}
/*** Fine funzione main ***/

/*** Inizion funzione per generare ed inizializzare la matrice ***/
int generate_matrix(char *matrix, int O_pct, int X_pct) {
    int row, column, random;

    srand(time(NULL) + MASTER);     // Genera un numero casuale

    // Controllo sulla grandezza della matrice
    if (ROWS <= 0 || COLUMNS <= 0) {
        printf("\033[1;31mERRORE\033[0m! La matrice deve essere più grande di 0 * 0\n\n");
        return 0;
    }

    // Controllo sulla percentuale di agenti rossi e blu
    if (O_pct + X_pct >= 100) {
        printf("\033[1;31mERRORE\033[0m! Rosso: %d%%, Blu: %d%%. La somma non può essere >= 100.\n\n", O_pct, X_pct);
        return 0;
    }

    for (row = 0; row < ROWS; row++) {
        for (column = 0; column < COLUMNS; column++) {
            random = rand() % 100;       // Numero casuale tra 0 e 99

            if ((random >= 0) && (random < O_pct)) {
                *(matrix + (row * COLUMNS) + column) = AGENT_O;
            }

            if ((random >= O_pct) && (random < O_pct + X_pct)) {
                *(matrix + (row * COLUMNS) + column) = AGENT_X;
            }

            if ((random >= O_pct + X_pct) && (random < 100)) {
                *(matrix + (row * COLUMNS) + column) = EMPTY;
            }
        }
    }

    return 1;
}
/*** Fine funzione per generare ed inizializzare la matrice ***/

/*** Inizio funzione per suddividere la matrice tra i processi ***/
int subdivide_matrix(int world_size, int *displacements, int *sendcounts, int *rows_per_process) {
    // Controllo sulla correttezza dell'input
    if (world_size <= 0 || ROWS <= 0 || COLUMNS <= 0 || displacements == NULL || sendcounts == NULL || rows_per_process == NULL) {
        printf("\033[1;31mERROR\033[0m! Invalid input in subdivide_matrix.\n\n");
        return 0;
    }

    int divisione = (ROWS) / (world_size);     // Numero di righe da asseganre ad ogni processo
    int resto = (ROWS) % (world_size);         // Nel caso in cui il numero di righe non sia divisibile per il nuermo di processori
    int displacement = 0;                      // All'inizio viene settato a 0

    // Assegna il giusto numero di righe ad ogni processo
    for (int i = 0; i < world_size; i++) {
        sendcounts[i] = divisione * COLUMNS;
        rows_per_process[i] = divisione;

        if (resto > 0) {
            sendcounts[i] += COLUMNS;  // Assegna una riga in più
            rows_per_process[i]++;
            resto--;
        }

        displacements[i] = displacement;
        displacement += sendcounts[i];

        // Nel calcolo vengono assegnate delle righe in più per contenere quelle dei processi adiacenti durante lo scambio
        rows_per_process[i] += (i == 0 || i == world_size - 1) ? 1 : 2;
    }

    return 1;
}
/*** Fine funzione per suddividere la matrice tra i processi ***/

/*** Inizio funzione per scambiare le righe dei processi vicini ***/
void exchange_rows(int rank, int world_size, int original_rows, char *sub_matrix, MPI_Comm communicator) {
    int neighbour_up, neighbour_down;
    MPI_Status status;            // Ritorna lo stato di una operazione di recive, wait o test
    MPI_Request request_up;       // Per capire quando un operazione non bloccate è finita
    MPI_Request request_down;     // Per capire quando un operazione non bloccate è finita

    // Rappresentano le righe dei processi adiacenti
    neighbour_up = (rank + 1) % world_size;
    neighbour_down = (rank + world_size - 1) % world_size;

    int my_last_row_pos = (original_rows - 1) * COLUMNS;                           // Indice della riga che l'ultiomo processo deve mandare al vicino superiore
    int neighbour_down_row_pos = original_rows * COLUMNS;                          // Indice della riga dove sarà salvata l'ultima riga del vicino precedente
    int neighbour_up_row_pos = (original_rows + ((rank == 0) ? 0 : 1)) * COLUMNS;  //Indice della riga dove sarà la prima riga del vicino superiore

    // Scambia le righe con i processi adiacenti
    if (rank != 0) {
        MPI_Isend(sub_matrix, COLUMNS, MPI_CHAR, neighbour_down, 99, communicator, &request_up);                                // (sendbuf, count, dtatype, destbuf, tag, comm, out request)
        MPI_Irecv(sub_matrix + neighbour_down_row_pos, COLUMNS, MPI_CHAR, neighbour_down, 99, communicator, &request_up);       // (recvbuff, count, datatype, sendbuff, tag, comm, out status)
    }

    // Scambia le righe con i processi adiacenti
    if (rank != world_size - 1) {
        MPI_Isend(sub_matrix + my_last_row_pos, COLUMNS, MPI_CHAR, neighbour_up, 99, communicator, &request_down);
        MPI_Irecv(sub_matrix + neighbour_up_row_pos, COLUMNS, MPI_CHAR, neighbour_up, 99, communicator, &request_down);
    }

    // Per completare la comunicazione non bloccante
    if (rank != 0) {
        MPI_Wait(&request_up, NULL);      // (request, status)
    }
    if (rank != world_size - 1) {
        MPI_Wait(&request_down, NULL);
    }
}
/*** Fine funzione per scambiare le righe dei processi vicini ***/

/*** Inizio funzione per calcolare gli agenti da spostare ***/
int *calculate_move(int rank, int world_size, int original_rows, int total_rows, char *sub_matrix, int *unsatisfied_agents) {
    // Alloca spazio per la matriche che conterra gli agenti che si vogliono sposare, inizialmente gli agenti insodisfatti sono 0 (ancora devono essere calcolati)
    int *mat = (int *)malloc(original_rows * COLUMNS * sizeof(int));
    *unsatisfied_agents = 0;

    // Cicla su tutta la matrice per calcolare gli agenti insodisfatti
    for (int i = 0; i < original_rows; i++) {
        for (int j = 0; j < COLUMNS; j++) {
            if (sub_matrix[i * COLUMNS + j] != EMPTY) {  // Se la cella non è vuota bisogna calcolare se lagente è sodisfatto o meno
                int satisfied = is_satisfied(rank, world_size, original_rows, total_rows, i * COLUMNS, j, sub_matrix);
                // is_satisfied = 1: non soddisfatto
                // is_satisfied = 0: soddisfatto
                mat[i * COLUMNS + j] = satisfied ? 0 : 1;

                // Se l'agente non è soddisfatto, viene incrementato il numero degli agenti insoddisfatti (servirà per il calcolo delle destinazioni possibili)
                if (mat[i * COLUMNS + j] == 1)
                    *unsatisfied_agents += 1;
            }
            // -1: la cella è vuota (è libera per chi vuole spostarsi)
            else
                mat[i * COLUMNS + j] = -1;
        }
    }

    return mat;
}
/*** Fine funzine per calcolare gli agenti che si vogliono spostare ***/

/*** Inizio funzione per calcolare se un agente è sodisfatto ***/
int is_satisfied(int rank, int world_size, int rows_size, int total_rows, int row, int column, char *sub_matrix) {
    int left_index, right_index;
    int ngh_precedent_row, ngh_next_row;  // Righe dei processi adiacenti
    char neighbours[8];                   // Matrice delle 8 celle adiacenti ad un agente

    char current_element = sub_matrix[row + column];
    int neighbours_count = 8;
    int similar = 0;

    // Calcolo dell'indice destro e sinistro (rappresentano le celle a destra e sinistra dell'agente)
    left_index = (COLUMNS + column - 1 != COLUMNS - 1) ? (COLUMNS + column - 1) % COLUMNS : -1;
    right_index = ((COLUMNS + column + 1) % COLUMNS != 0) ? (COLUMNS + column + 1) % COLUMNS : -1;

    // Calcolo la posizione delle righe precedenti e successive (rappresentano le celle sopra e sotto l'agente)
    if (rank == 0) {
        ngh_precedent_row = -1;                                // Non ci sono righe precedenti
        ngh_next_row = total_rows * COLUMNS - COLUMNS;         // Posizione della riga successiva all'ultima del processo
    } else if (rank == world_size - 1) {                       // Non ha righe successive
        ngh_precedent_row = total_rows * COLUMNS - COLUMNS;
        ngh_next_row = -1;                                     // Non ha righe successive
    } else {                                                   // Processi con righe centrali
        ngh_precedent_row = total_rows * COLUMNS - COLUMNS - COLUMNS;
        ngh_next_row = total_rows * COLUMNS - COLUMNS;
    }

    if (row != 0) {
        if (left_index != -1)                                            // L'elemento a sinistra esiste
            neighbours[0] = sub_matrix[row - COLUMNS + left_index];
        else                                                             // L'elemento a sinistra non esiste
            neighbours[0] = '\0';
    } else {
        if (left_index != -1) {
            neighbours[0] = rank == 0 ? '\0' : sub_matrix[ngh_precedent_row + left_index];
        } else
            neighbours[0] = '\0';
    }

    if (row != 0) {
        neighbours[1] = sub_matrix[row - COLUMNS + column];
    } else {
        neighbours[1] = rank == 0 ? '\0' : sub_matrix[ngh_precedent_row + column];
    }

    if (row != 0) {
        if (right_index != -1)                                          // L'elemento a destra esiste
            neighbours[2] = sub_matrix[row - COLUMNS + right_index];
        else                                                            // L'elemento a destra non esiste
            neighbours[2] = '\0';
    } else {
        if (right_index != -1)
            neighbours[2] = rank == 0 ? '\0' : sub_matrix[ngh_precedent_row + right_index];
        else
            neighbours[2] = '\0';
    }

    // Riga corrente
    neighbours[3] = left_index != -1 ? sub_matrix[row + left_index] : '\0';
    neighbours[4] = right_index != -1 ? sub_matrix[row + right_index] : '\0';

    // Riga successiva
    if (row != rows_size - 1) {
        if (left_index != -1)                                          // L'elemento a sinistra esiste
            neighbours[5] = sub_matrix[row + COLUMNS + left_index];
        else                                                           // L'elemento a sinistra non esiste
            neighbours[5] = '\0';
    } else {
        if (left_index != -1) {
            neighbours[5] = rank == world_size - 1 ? '\0' : sub_matrix[ngh_next_row + left_index];
        } else
            neighbours[5] = '\0';
    }

    if (row != rows_size - 1) {
        neighbours[6] = sub_matrix[row + COLUMNS + column];
    } else {
        neighbours[6] = rank == world_size - 1 ? '\0' : sub_matrix[ngh_next_row + column];
    }

    if (row != rows_size - 1) {
        if (right_index != -1)                                         // L'elemento a destra esiste
            neighbours[7] = sub_matrix[row + COLUMNS + right_index];
        else                                                           // L'elemento a destra non esiste
            neighbours[7] = '\0';
    } else {
        if (right_index != -1)
            neighbours[7] = rank == world_size - 1 ? '\0' : sub_matrix[ngh_next_row + right_index];
        else
            neighbours[7] = '\0';
    }

    // Conta degli elementi adiacenti all'agente simili ad esso
    for (int i = 0; i < 8; i++) {
        if (neighbours[i] == current_element)
            similar++;
        else if (neighbours[i] == '\0')
            neighbours_count--;
    }

    // Calcolo della soddisfazione
    if ((((double)100 / neighbours_count) * similar) >= SAT_PERCENTAGE)
        return 1;
    else
        return 0;
}
/*** Fine funzione per calcolare se un agente è sodisfatto **/

/*** Inizio funzione per calcolare il numero di celle vuote locali ad un processo ***/
voidCell *calculate_local_void_cells(int original_rows, char *sub_matrix, int displacement, int *local_void_cells) {
    voidCell *void_cells;    // Array che contiene le celle vuote ([riga][colonna]) -> ([riga * COLUMNS + colonna])
    int ind = 0;             // Indice che indica quante celle vuote ha trovato. Alla fine lo assegna a 'local_void_cells'

    void_cells = malloc(original_rows * COLUMNS * sizeof(voidCell));

    // Calcolo delle celle vuote
    for (int i = 0; i < original_rows; i++)
        for (int j = 0; j < COLUMNS; j++)
            if (sub_matrix[i * COLUMNS + j] == EMPTY) {
                voidCell temp = {(displacement + i * COLUMNS), j};
                void_cells[ind] = temp;
                ind++;
            }

    void_cells = realloc(void_cells, ind * sizeof(voidCell));
    *local_void_cells = ind;

    return void_cells;
}
/*** Fine funzione per calcolare il numero di celle vuote locali ad un processo ***/

/*** Inizio funzione per unire tutte le celle vuote dei processi e restituire quelle di destinazione per il processo i-esimo ***/
voidCell *assign_void_cells(int rank, int world_size, int number_of_local_void_cells, voidCell *local_void_cells, int *number_of_void_cells_to_return, MPI_Datatype datatype, int unsatisfied_agents) {
    int number_of_global_void_cells[world_size];     // Array che contiene il numero di celle vuote per ogni processo
    int displacements[world_size];                   // Displacements per la prima gather (displs[rank] == numero di celle vuote del rank)
    int number_of_total_void_cells = 0;              // Numero di celle vuote in tutta la matrice
    voidCell *global_void_cells;                     // Array con tutte le posizioni delle celle vuote - es. { [0][1], [3][2], [5][8] }
    int *void_cells_per_process;                     // Array che indica quante celle vuote vengono assegnate ad un processo
    int global_unsatisfied_agents[world_size];       // Array che contiene il numero degli agenti insoddisfatti per ogni processo

    global_void_cells = malloc(ROWS * COLUMNS * sizeof(voidCell));
    void_cells_per_process = malloc(ROWS * COLUMNS * sizeof(voidCell));

    // Il numero di celle vuote di ogni processo viene condiviso con tutti gli altri
    MPI_Allgather(&number_of_local_void_cells, 1, MPI_INT, number_of_global_void_cells, 1, MPI_INT, MPI_COMM_WORLD);       // (sendbuff, sendcount, datatype, destbuff, destcount, datatype, comm)

    // Calcolo del displacement e del numero totale di celle vuote
    for (int i = 0; i < world_size; i++) {
        displacements[i] = i == 0 ? 0 : displacements[i - 1] + number_of_global_void_cells[i - 1];
        number_of_total_void_cells += number_of_global_void_cells[i];
    }

    // Vengono raggruppate tutte le celle vuote
    MPI_Allgatherv(local_void_cells, number_of_local_void_cells, datatype, global_void_cells, number_of_global_void_cells, displacements, datatype, MPI_COMM_WORLD);        // (sendbuff, sendcount, senddatatype, destbuff, destcount, displacements, destdatatype, comm)

    // Il numero di agenti insoddisfatti per ogni processo viene messo in un array
    MPI_Allgather(&unsatisfied_agents, 1, MPI_INT, global_unsatisfied_agents, 1, MPI_INT, MPI_COMM_WORLD);     // (sendubb, sendcount, datatype, destbuff, destcount, datatype, comm)

    // L'array con le celle vuote viene mescolato, viene usato lo stesso seme per ogni processo
    srand(SEED);                      // Inizzializza il seme
    for (int i = 0; i < number_of_total_void_cells; i++) {
        int destination = ((rand() % (number_of_total_void_cells)) + 0);
        voidCell tmp = global_void_cells[destination];
        global_void_cells[destination] = global_void_cells[i];
        global_void_cells[i] = tmp;
    }

    // Calcolo le posizioni locali al processo e Vengono assegnate le celle vuote ai processi
    int divisione = number_of_total_void_cells / world_size;           // Celle vuote da asseganre ad ogni processo
    int resto = number_of_total_void_cells % world_size;               // Resto se != 0 bisogna assegnare più celle vuote ad un processo
    int displacement = 0;
    for (int i = 0; i < world_size; i++) {
        void_cells_per_process[i] = divisione > global_unsatisfied_agents[i] ? global_unsatisfied_agents[i] : divisione;

        if (resto > 0 && !(divisione > global_unsatisfied_agents[i])) {
            void_cells_per_process[i]++;
            resto--;
        }

        displacements[i] = displacement;
        displacement += void_cells_per_process[i];
    }

    // Ad ogni processo viene assegnato un numero di celle vuote
    *number_of_void_cells_to_return = void_cells_per_process[rank];
    voidCell *toReturn = malloc(sizeof(voidCell) * void_cells_per_process[rank]);      // Contiene le celle vuote da assegnare ad ogni processo
    MPI_Scatterv(global_void_cells, void_cells_per_process, displacements, datatype, toReturn, void_cells_per_process[rank], datatype, MASTER, MPI_COMM_WORLD);    // (sendbuff, sendcount, displacements, datatype, destbuff, destcount, datatype, root, comm)

    free(global_void_cells);
    free(void_cells_per_process);

    return toReturn;
}
/*** Fine funzione per unire tutte le celle vuote dei processi e restituire quelle di destinazione per il processo i-esimo ***/

/*** Inizo funzione per calcolare a quale processo appartiene una determinata riga della matrice ***/
int calculate_source(int world_size, int *displacement, int *sendcounts, int row) {
    int toReturn = 0;

    // Si calcola a quale processo (rank) appartiene una cella di destinazione
    for (int tmp_rank = 0; tmp_rank < world_size; tmp_rank++) {
        int start = displacement[tmp_rank] / COLUMNS;
        int finish = (start + sendcounts[tmp_rank] / COLUMNS);

        if (row >= start && row < finish) {
            toReturn = tmp_rank;
            break;
        }
    }
    return toReturn;
}
/*** Fine funzione per calcolare a quale processo appartiene una determinata roga della matrice ***/

/*** Inizio funzione per spostare gli agenti ***/
void move(int rank, int world_size, int original_rows, char *sub_matrix, int *want_move, voidCell *destinations, int num_assigned_void_cells, int *displacements, int *sendcounts, MPI_Datatype move_agent_type) {
    int num_elems_to_send_to[world_size];      // Array che contiene il numero di moveAgent da mandare al processo i-esimo
    int used_void_cells_assigned = 0;          // Il numero delle celle vuote che sono state assegnate al processo e che ha usato.
    moveAgent **data;                          // Matrice che contiene sulle righe i processi e sulle colonne la cella di destinazione dell'agente che vuole spostarsi

    memset(num_elems_to_send_to, 0, sizeof(num_elems_to_send_to));                    // Setta a 0 gli elementi dell'array appena creato
    data = (moveAgent **)malloc(sizeof(moveAgent *) * world_size);                    // Alloca spazio per data
    for (int i = 0; i < world_size; i++)
        data[i] = (moveAgent *)malloc(sizeof(moveAgent) * num_assigned_void_cells);   // Alloca spazio per ogni cella

    // Si itera finche non finiscono le righe a disposizione o il numero di celle vuote
    for (int i = 0; i < original_rows && used_void_cells_assigned < num_assigned_void_cells; i++) {
        for (int j = 0; j < COLUMNS && used_void_cells_assigned < num_assigned_void_cells; j++) {
            // Sposta l'agente in una cella libera
            if (want_move[i * COLUMNS + j] == 1) {                                                                          // Se l'agente vuole spostarsi
                voidCell destination = destinations[used_void_cells_assigned];                                              // Gli viene assegnata una cella libera
                int receiver = calculate_source(world_size, displacements, sendcounts, destination.row_index / COLUMNS);    // Si verifica a che processo appartiene la cella di destinazionr

                // La cella di destinazione appartiene al processo stesso, l'agente viene subito spostato
                if (receiver == rank) {
                    int startRow = displacements[rank];                                                 // Riga iniziale
                    int destRow = destination.row_index - startRow;                                     // Riga di destinazione

                    sub_matrix[destRow + destination.column_index] = sub_matrix[i * COLUMNS + j];       // Sposta l'agente
                    sub_matrix[i * COLUMNS + j] = EMPTY;                                                // Libera lo spazio nella sottomatrice

                    want_move[destRow + destination.column_index] = 0;                                  // Non rendere più disponibile lo spazio disponibile per altri
                    want_move[i * COLUMNS + j] = -1;                                                    // Libera questo spazio precedente
                }
                // La cella di destinazione non appartiene al processo stesso
                else {
                    int startRow = displacements[receiver];                                             // Riga iniziale del destinatario
                    int destRow = destination.row_index - startRow;                                     // Riga di destinazione del destinatario

                    moveAgent var = {destRow, destination.column_index, sub_matrix[i * COLUMNS + j]};   // Informazioni e agente che vuole spostarsi
                    data[receiver][num_elems_to_send_to[receiver]] = var;                               // Setta al processo 'receiver' la X-esima colonna con la cella di destinazione dell'agente
                    num_elems_to_send_to[receiver] += 1;                                                // Aggiorna il numero di elementi che deve mandare al processo 'receiver'

                    sub_matrix[i * COLUMNS + j] = EMPTY;                                                // Libera lo spazio nella sottomatrice
                    want_move[i * COLUMNS + j] = -1;                                                    // Libera questo spazio precedente
                }

                used_void_cells_assigned++;                                                             // Aggiorna il numero di celle vuote che ha usato
            }
        }
    }

    // Tutti i processi vengono sincronizzati
    synchronize(rank, world_size, num_elems_to_send_to, num_assigned_void_cells, data, original_rows, sub_matrix, move_agent_type);
}
/*** Fine funzione per spostare gli agenti ***/

/*** Inizio funzione per sincronizzare gli postamenti tra i processi ***/
void synchronize(int rank, int world_size, int *num_elems_to_send_to, int num_assigned_void_cells, moveAgent **data, int original_rows, char *sub_matrix, MPI_Datatype move_agent_type) {
    int my_void_cell_used_by[world_size];    // Array che contiene in ogni cella il numero di elementi che il processo i-esimo vuole scrivere nelle celle della sottomatrice
    MPI_Request requests1[world_size];       // Array per le prime MPI_Irecv e MPI_Wait
    MPI_Request requests2[world_size];       // Array per le seconde MPI_Irecv e MPI_Wait
    moveAgent **moved_agents;                // Matrice degli agenti che il processo ha ricevuto che deve aggiornare nella sottomatrice
    moveAgent **elements_to_send;            // Array che contiene gli elementi da mandare al processo i-esimo

    moved_agents = (moveAgent **)malloc(sizeof(moveAgent *) * world_size - 1);
    elements_to_send = (moveAgent **)malloc(sizeof(moveAgent *) * world_size - 1);

    // Vengono calcolate quante celle sono state usate dei num_elems_to_send_to
    for (int i = 0; i < world_size; i++) {
        if (i == rank) continue;

        MPI_Isend(&num_elems_to_send_to[i], 1, MPI_INT, i, 99, MPI_COMM_WORLD, &requests1[i]);  // Manda al processo i il numero di celle che ha usato
        MPI_Irecv(&my_void_cell_used_by[i], 1, MPI_INT, i, 99, MPI_COMM_WORLD, &requests1[i]);  // Riceve dal processo i il numero di celle del processo che lui ha usato
    }

    // Manda/riceve al/dal processo i-esimo tutte le celle di destinazione dove deve scrivere/salvare i suoi agenti
    for (int i = 0; i < world_size; i++) {
        if (i == rank) continue;

        int number_of_elems_to_send;  // Numero di elementi da mandare al processo i-esimo
        number_of_elems_to_send = num_elems_to_send_to[i];

        elements_to_send[i] = (moveAgent *)malloc(number_of_elems_to_send * sizeof(moveAgent));  // Alloca spazio per la i-esima riga

        for (int j = 0; j < number_of_elems_to_send; j++)
            elements_to_send[i][j] = data[i][j];

        MPI_Wait(&requests1[i], NULL);  // Aspetta che la prima Irecv riceva il numero di elementi che gli altri processi vogliono scrivere nelle sue celle della sottomatrice

        moved_agents[i] = (moveAgent *)malloc(my_void_cell_used_by[i] * sizeof(moveAgent));
        MPI_Isend(elements_to_send[i], number_of_elems_to_send, move_agent_type, i, 100, MPI_COMM_WORLD, &requests2[i]);
        MPI_Irecv(moved_agents[i], my_void_cell_used_by[i], move_agent_type, i, 100, MPI_COMM_WORLD, &requests2[i]);
    }

    // Aspetta che le seconde MPI_Wait terminino
    for (int i = 0; i < world_size; i++) {
        if (i == rank) continue;

        MPI_Wait(&requests2[i], NULL);
    }

    // Scrive gli agenti 'nuovi' nelle celle di destinazione
    for (int i = 0; i < world_size; i++) {
        if (i == rank) continue;

        for (int k = 0; k < my_void_cell_used_by[i]; k++)
            sub_matrix[moved_agents[i][k].destination_row + moved_agents[i][k].destination_column] = moved_agents[i][k].agent;  // Scrive l'agente nella cella vuota
    }

    // Dealloca
    for (int i = 0; i < world_size; i++) {
        free(data[i]);
        if (i == rank) continue;

        free(moved_agents[i]);
    }
    free(data);
    free(moved_agents);
    free(elements_to_send);
}
/*** Fine funzioe per sincronizzare gli spostamenti tra i processi ***/

/*** Inizio funzione per definire il tipo voidCell ***/
void define_voidCell_type(MPI_Datatype *VOID_CELL_TYPE) {
    int vc_items = 2;                   // Numero di items
    int vc_block_length[2] = {1, 1};    // Lunghezza dei blocchi

    MPI_Aint vc_offsets[2];    // Definisce il tipo di una variabile che contiene un indirizzo di memoria
    voidCell void_cell;        // Di tipo void_cell
    MPI_Aint vc_base_address;
    MPI_Get_address(&void_cell, &vc_base_address);      // Restituisce l'indirizzo di una locazione in memoria (location, addres)
    MPI_Get_address(&void_cell.row_index, &vc_offsets[0]);
    MPI_Get_address(&void_cell.column_index, &vc_offsets[1]);
    vc_offsets[0] = MPI_Aint_diff(vc_offsets[0], vc_base_address);   // Calcola l'offset fscendo la differenza tra due indirizzi
    vc_offsets[1] = MPI_Aint_diff(vc_offsets[1], vc_base_address);

    MPI_Datatype vc_types[2] = {MPI_INT, MPI_INT};
    MPI_Type_create_struct(vc_items, vc_block_length, vc_offsets, vc_types, VOID_CELL_TYPE);   // Crea un MPI datatype (numero di blocchi da creare, array che contiene la lunghezza di ogni blocco, displacment per ogni blocco, array che contiene l'MPI datatype da replicare per ogni blocco, variabile in cui mettere il nuovo datatype)
    MPI_Type_commit(VOID_CELL_TYPE);            // Commit del tipo
}
/*** Fine funzie per definire il tipo voidCell ***/

/*** Inizio funzione per definire il tipo moveAgent ***/
void define_moveAgent_type(MPI_Datatype *MOVE_AGENT_TYPE) {
    int ma_items = 3;                    // Numero di items
    int ma_block_length[3] = {1, 1, 1};  // Lunghezza dei blocchi

    MPI_Aint ma_offsets[3];      // Definisce il tipo di una variabile che contiene un indirizzo di memoria
    moveAgent move_agent;        // Di tipo move_agent
    MPI_Aint ma_base_address;
    MPI_Get_address(&move_agent, &ma_base_address);     // Restituisce l'indirizzo di una locazione in memoria (location, addres)
    MPI_Get_address(&move_agent.destination_row, &ma_offsets[0]);
    MPI_Get_address(&move_agent.destination_column, &ma_offsets[1]);
    MPI_Get_address(&move_agent.agent, &ma_offsets[2]);
    ma_offsets[0] = MPI_Aint_diff(ma_offsets[0], ma_base_address);    // Calcola l'offset fscendo la differenza tra due indirizzi
    ma_offsets[1] = MPI_Aint_diff(ma_offsets[1], ma_base_address);
    ma_offsets[2] = MPI_Aint_diff(ma_offsets[2], ma_base_address);

    MPI_Datatype ma_types[3] = {MPI_INT, MPI_INT, MPI_CHAR};
    MPI_Type_create_struct(ma_items, ma_block_length, ma_offsets, ma_types, MOVE_AGENT_TYPE);    // Crea un MPI datatype (numero di blocchi da creare, array che contiene la lunghezza di ogni blocco, displacment per ogni blocco, array che contiene l'MPI datatype da replicare per ogni blocco, variabile in cui mettere il nuovo datatype)
    MPI_Type_commit(MOVE_AGENT_TYPE);                       // Commit del tipo
}
/*** Fine funzioe per definire il tipo moveAgent ***/

/*** Inizio funzione per calcolare la soddisfazione finale di tutti gli agenti ***/
void calculate_total_satisfaction(int rank, int world_size, char *matrix) {
    int total_agents = 0;                                                                 // Agenti totali
    int satisfied_agents = 0;                                                             // Agenti non soddisfatti
    moveAgent *not_satisfied_agents = malloc(ROWS * COLUMNS * sizeof(moveAgent));         // Agenti soddisfatti

    int index = 0;

    for (int i = 0; i < ROWS; i++)
        for (int j = 0; j < COLUMNS; j++)
            if (matrix[i * COLUMNS + j] != EMPTY) {
                total_agents++;
                if (is_satisfied(rank, world_size, ROWS, ROWS, i * COLUMNS, j, matrix)) {
                    satisfied_agents++;
                } else {
                    moveAgent var = {i * COLUMNS, j, matrix[i * COLUMNS + j]};
                    not_satisfied_agents[index] = var;
                    index++;
                }
            }

    // Rapporto tra agenti soddisfatti e insodisfatti
    float average = ((double)satisfied_agents / (double)total_agents) * 100;
    printf("\nInfo:\n");
    printf("- Agenti totali: %d\n", total_agents);
    printf("- Agenti soddisfatti: %d\n", satisfied_agents);
    if (index != 0) {
        printf("- Agenti non soddisfatti: %d\n", index);
    }
    printf("Percentuale di soddisfazione: %.3f%%\n", average);

    free(not_satisfied_agents);
}
/*** Fine funzione per calcolare la soddisfazione finale di tutti gli agenti ***/

/*** Inizio funzione per visualizzare la matrice ***/
void print_matrix(int rows_size, int column_size, char *matrix) {
    int i;
    if (rows_size <= 0 || column_size <= 0 || matrix == NULL) {
        printf("\033[1;31mERROR\033[0m! Invalid input.\n");
        MPI_Abort(MPI_COMM_WORLD, 0);
        return;
    }

    for (i = 0; i < rows_size * column_size; i++) {
        if (*(matrix + i) == AGENT_X) {
            printf(BLUE("%c") " ", *(matrix + i));
        } else if (*(matrix + i) == AGENT_O) {
            printf(RED("%c") " ", *(matrix + i));
        } else
            printf("%c ", *(matrix + i));

        if ((i + 1) % column_size == 0 && i != 0)
            printf("\n");
    }
}
/*** Fine funzione per visualizzare la matrice ***/

/*** Inizio funzione per terminare in caso di errori ***/
void err_finish(int *sendcounts, int *displacements, int *rows_per_process) {
    free(sendcounts);
    free(displacements);
    free(rows_per_process);
    MPI_Abort(MPI_COMM_WORLD, MPI_ERR_COUNT);
    exit(0);
}
/*** Fine funzione per terminare in caso di errori ***/

/*** Inizio funzioe di demo con matrice statica ***/
void test_init_matrix(char *matrix, int O_pct, int X_pct) {
    int row, column, random;

    matrix[0 * COLUMNS + 0] = AGENT_X;
    matrix[0 * COLUMNS + 1] = AGENT_O;
    matrix[0 * COLUMNS + 2] = AGENT_X;
    matrix[0 * COLUMNS + 3] = AGENT_X;
    matrix[0 * COLUMNS + 4] = AGENT_X;
    matrix[0 * COLUMNS + 5] = AGENT_O;
    matrix[0 * COLUMNS + 6] = AGENT_O;
    matrix[0 * COLUMNS + 7] = EMPTY;
    matrix[0 * COLUMNS + 8] = AGENT_O;
    matrix[0 * COLUMNS + 9] = AGENT_X;
    matrix[1 * COLUMNS + 0] = AGENT_X;
    matrix[1 * COLUMNS + 1] = AGENT_X;
    matrix[1 * COLUMNS + 2] = AGENT_X;
    matrix[1 * COLUMNS + 3] = AGENT_O;
    matrix[1 * COLUMNS + 4] = AGENT_O;
    matrix[1 * COLUMNS + 5] = EMPTY;
    matrix[1 * COLUMNS + 6] = EMPTY;
    matrix[1 * COLUMNS + 7] = AGENT_X;
    matrix[1 * COLUMNS + 8] = AGENT_X;
    matrix[1 * COLUMNS + 9] = EMPTY;
    matrix[2 * COLUMNS + 0] = AGENT_X;
    matrix[2 * COLUMNS + 1] = AGENT_O;
    matrix[2 * COLUMNS + 2] = EMPTY;
    matrix[2 * COLUMNS + 3] = AGENT_O;
    matrix[2 * COLUMNS + 4] = AGENT_O;
    matrix[2 * COLUMNS + 5] = AGENT_X;
    matrix[2 * COLUMNS + 6] = AGENT_O;
    matrix[2 * COLUMNS + 7] = AGENT_X;
    matrix[2 * COLUMNS + 8] = AGENT_O;
    matrix[2 * COLUMNS + 9] = AGENT_X;
    matrix[3 * COLUMNS + 0] = AGENT_O;
    matrix[3 * COLUMNS + 1] = EMPTY;
    matrix[3 * COLUMNS + 2] = EMPTY;
    matrix[3 * COLUMNS + 3] = AGENT_O;
    matrix[3 * COLUMNS + 4] = AGENT_O;
    matrix[3 * COLUMNS + 5] = EMPTY;
    matrix[3 * COLUMNS + 6] = AGENT_X;
    matrix[3 * COLUMNS + 7] = EMPTY;
    matrix[3 * COLUMNS + 8] = AGENT_O;
    matrix[3 * COLUMNS + 9] = EMPTY;
    matrix[4 * COLUMNS + 0] = EMPTY;
    matrix[4 * COLUMNS + 1] = AGENT_O;
    matrix[4 * COLUMNS + 2] = AGENT_O;
    matrix[4 * COLUMNS + 3] = AGENT_X;
    matrix[4 * COLUMNS + 4] = AGENT_O;
    matrix[4 * COLUMNS + 5] = AGENT_X;
    matrix[4 * COLUMNS + 6] = EMPTY;
    matrix[4 * COLUMNS + 7] = AGENT_X;
    matrix[4 * COLUMNS + 8] = AGENT_O;
    matrix[4 * COLUMNS + 9] = EMPTY;
    matrix[5 * COLUMNS + 0] = AGENT_X;
    matrix[5 * COLUMNS + 1] = AGENT_O;
    matrix[5 * COLUMNS + 2] = AGENT_X;
    matrix[5 * COLUMNS + 3] = AGENT_X;
    matrix[5 * COLUMNS + 4] = AGENT_X;
    matrix[5 * COLUMNS + 5] = AGENT_O;
    matrix[5 * COLUMNS + 6] = EMPTY;
    matrix[5 * COLUMNS + 7] = AGENT_X;
    matrix[5 * COLUMNS + 8] = AGENT_O;
    matrix[5 * COLUMNS + 9] = EMPTY;
    matrix[6 * COLUMNS + 0] = AGENT_X;
    matrix[6 * COLUMNS + 1] = AGENT_X;
    matrix[6 * COLUMNS + 2] = AGENT_X;
    matrix[6 * COLUMNS + 3] = AGENT_O;
    matrix[6 * COLUMNS + 4] = AGENT_O;
    matrix[6 * COLUMNS + 5] = AGENT_X;
    matrix[6 * COLUMNS + 6] = AGENT_X;
    matrix[6 * COLUMNS + 7] = EMPTY;
    matrix[6 * COLUMNS + 8] = AGENT_O;
    matrix[6 * COLUMNS + 9] = EMPTY;
    matrix[7 * COLUMNS + 0] = AGENT_X;
    matrix[7 * COLUMNS + 1] = AGENT_O;
    matrix[7 * COLUMNS + 2] = EMPTY;
    matrix[7 * COLUMNS + 3] = AGENT_O;
    matrix[7 * COLUMNS + 4] = AGENT_O;
    matrix[7 * COLUMNS + 5] = AGENT_O;
    matrix[7 * COLUMNS + 6] = AGENT_X;
    matrix[7 * COLUMNS + 7] = AGENT_O;
    matrix[7 * COLUMNS + 8] = EMPTY;
    matrix[7 * COLUMNS + 9] = AGENT_X;
    matrix[8 * COLUMNS + 0] = AGENT_O;
    matrix[8 * COLUMNS + 1] = EMPTY;
    matrix[8 * COLUMNS + 2] = EMPTY;
    matrix[8 * COLUMNS + 3] = AGENT_O;
    matrix[8 * COLUMNS + 4] = AGENT_O;
    matrix[8 * COLUMNS + 5] = EMPTY;
    matrix[8 * COLUMNS + 6] = AGENT_X;
    matrix[8 * COLUMNS + 7] = AGENT_O;
    matrix[8 * COLUMNS + 8] = AGENT_X;
    matrix[8 * COLUMNS + 9] = EMPTY;
    matrix[9 * COLUMNS + 0] = EMPTY;
    matrix[9 * COLUMNS + 1] = AGENT_O;
    matrix[9 * COLUMNS + 2] = AGENT_O;
    matrix[9 * COLUMNS + 3] = AGENT_X;
    matrix[9 * COLUMNS + 4] = AGENT_O;
    matrix[9 * COLUMNS + 5] = AGENT_X;
    matrix[9 * COLUMNS + 6] = AGENT_X;
    matrix[9 * COLUMNS + 7] = AGENT_O;
    matrix[9 * COLUMNS + 8] = AGENT_O;
    matrix[9 * COLUMNS + 9] = AGENT_X;
}
/*** Fine funzione di demo con matroce statica ***/