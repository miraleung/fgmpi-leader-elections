/**
 * lcr-passthru.c
 * 
 * @author Mira Leung 
 *
 * March 28, 2014
 *
 * An implementation of Lelann/Chang-Roberts', except that it checks for the 
 * minimum uid seen so far, instead of against its own.
 * Unidirectional ring.
 *
 * uids are not randomly assigned, so that we can have a single initiator.
 *
 */

#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fgmpi.h>

// Tags
#define TAG_PHASE1 2
#define TAG_ELECTION 3
#define TAG_NRECV 4
#define TAG_NSENT 5
#define TAG_MSGNUM 6

#define SIZE_MSG 2

// Process states
typedef enum { NONACTIVE, ACTIVE, LEADER } process_state; // A NONACTIVE process lost the election

int ceiling_log2(unsigned long long x);
int gcd(int size, int pnum);

/** FG-MPI Boilerplate begins **/
int lcr_passthru(int argc, char* argv[]);
FG_ProcessPtr_t binding_func(int argc, char** argv, int rank) {
  return (&lcr_passthru);
}

FG_MapPtr_t map_lookup(int argc, char** argv, char* str) {
  return (&binding_func);
}

int main(int argc, char *argv[]) {
  FGmpiexec(&argc, &argv, &map_lookup);
  return 0;
}

/** FG-MPI Boilerplate ends **/



/**
 * Main
 */
int lcr_passthru(int argc, char *argv[]) {

  if (argc != 2 && argc != 3) {
    printf("Usage: ./lcr-passthru [ -v ] <Process number>\n");
    exit(1);
  }

  int rank, size, uid, pnum;
  int tag, max_so_far;
  int recv_buf[2];

  int lnum_sent = 0, lnum_recv = 0;
  int tnum_sent = 0, tnum_recv = 0;

  process_state my_state = ACTIVE;

  int verbose = 0;
  if (argc == 3) {
    if (!strcmp(argv[1], "-v")) pnum = atoi(argv[2]), verbose = 1;
    else if (!strcmp(argv[2], "-v")) pnum = atoi(argv[1]),  verbose = 1;
  } else if (argc == 2) 
    pnum = atoi(argv[1]);
 

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Request request;
  MPI_Status status;

  
  if (pnum <= size || (int) pnum/size < 6 || gcd(size, pnum) != 1) {
    printf("Usage: pnum must be at least 6 times larger than and relatively coprime to size.\n");
    exit(1);
  }

 
  int send_neighbour = (rank+1) % size, recv_neighbour = rank - 1;
  if (!rank) recv_neighbour = size - 1;

  srand(time(NULL) + rank);
  //uid = (rand() % pnum);
  uid = ((rank+1)*pnum) % size;
  int initiator = (uid % size) == (size - 1)/2; 
  int participant = 0;
  int rnd = rand() % size;
  int canParticipate = (rnd % 5);


  max_so_far = uid;
  tag = TAG_PHASE1;

  if (initiator) {
    printf("Process %d is an initiator\n", rank);
    participant = 1;
    canParticipate = 1;
    MPI_Isend(&max_so_far, SIZE_MSG, MPI_INT, send_neighbour, tag, MPI_COMM_WORLD, &request);
    lnum_sent++;
  }

  //  Everyone is an initiator by default
  while (my_state == ACTIVE && canParticipate) { 
      MPI_Recv(&recv_buf, SIZE_MSG, MPI_INT, recv_neighbour, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      lnum_recv++;
 
      // Got an election message or a smaller uid than the least seen so far, so I know I lost
      if (status.MPI_TAG == TAG_ELECTION || recv_buf[0] > max_so_far) {
        max_so_far = recv_buf[0];
        my_state = NONACTIVE; // lost the election
        // forward the message, and break;
        MPI_Isend(recv_buf, SIZE_MSG, MPI_INT, send_neighbour, status.MPI_TAG, MPI_COMM_WORLD, &request);
        lnum_sent++;
        break; 
      }

      // Got my own uid back; I'm the leader
      if (recv_buf[0] == uid) {
        max_so_far = uid;
        my_state = LEADER;
        tag = TAG_ELECTION;
        MPI_Isend(&uid, SIZE_MSG, MPI_INT, send_neighbour, tag, MPI_COMM_WORLD, &request);
        lnum_sent++;
      } else if (recv_buf[0] < uid && !participant) {
        participant = 1;
        MPI_Isend(&uid, SIZE_MSG, MPI_INT, send_neighbour, TAG_PHASE1, MPI_COMM_WORLD, &request);
        lnum_sent++;
      }
  }

  // Non-candidates forward messages
  while (1) {
    MPI_Recv(recv_buf, SIZE_MSG, MPI_INT, recv_neighbour, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    lnum_recv++;
    if ((my_state == NONACTIVE || !canParticipate) && status.MPI_TAG == TAG_ELECTION) {
      if (recv_buf[0] >  max_so_far) max_so_far = recv_buf[0];
      MPI_Isend(recv_buf, SIZE_MSG, MPI_INT, send_neighbour, status.MPI_TAG, MPI_COMM_WORLD, &request);
      lnum_sent++;
      break;
    } else if  (my_state == LEADER && recv_buf[0] == uid && status.MPI_TAG == TAG_ELECTION) {
      if (recv_buf[0] > max_so_far) max_so_far = recv_buf[0];
      break;
    }

    MPI_Isend(recv_buf, SIZE_MSG, MPI_INT, send_neighbour, status.MPI_TAG, MPI_COMM_WORLD, &request);
    lnum_sent++;
  }

  if (canParticipate && participant && verbose) 
  printf("rank=%d, id=%d, leader=%d, mrcvd=%d, msent=%d\n", rank, uid, max_so_far == uid, lnum_recv, lnum_sent);

  // Non-leaders, send your local message totals
  int msgBuf[2] = { lnum_recv, lnum_sent };
  if (!canParticipate || !participant) msgBuf[0] = 0, msgBuf[1] = 0;
  if (my_state != LEADER) {
    MPI_Recv(recv_buf, SIZE_MSG, MPI_INT, recv_neighbour, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    if (status.MPI_TAG == TAG_MSGNUM && canParticipate && participant) msgBuf[0] += recv_buf[0], msgBuf[1] += recv_buf[1]; 
    else msgBuf[0] = recv_buf[0], msgBuf[1] = recv_buf[1];
  }
    MPI_Isend(msgBuf, SIZE_MSG, MPI_INT, send_neighbour, TAG_MSGNUM, MPI_COMM_WORLD, &request);

  // Leader receives/prints total number of messages sent and received
  if (my_state == LEADER && participant) {
      MPI_Recv(recv_buf, SIZE_MSG, MPI_INT, recv_neighbour, TAG_MSGNUM, MPI_COMM_WORLD, &status);      
      tnum_recv = recv_buf[0]; tnum_sent = recv_buf[1];

    printf("Leader: rank=%d, id=%d, trcvd=%d, tsent=%d\n", rank, uid, tnum_recv, tnum_sent);  
  }

  

  MPI_Finalize();
  return 0;
}


int gcd(int size, int pnum) {
  int k = size, m = pnum;  
  while (k != m) {
    if (k > m) k = k-m; 
    else m = m-k; 
   }
   return k;
}

int ceiling_log2(unsigned long long x) {
  static const unsigned long long t[6] = {
    0xFFFFFFFF00000000ull,
    0x00000000FFFF0000ull,
    0x000000000000FF00ull,
    0x00000000000000F0ull,
    0x000000000000000Cull,
    0x0000000000000002ull
  };

  int y = (((x & (x - 1)) == 0) ? 0 : 1);
  int j = 32, i;

  for (i = 0; i < 6; i++) {
    int k = (((x & t[i]) == 0) ? 0 : j);
    y += k, x >>= k, j >>= 1;
  }

  return y;
}

