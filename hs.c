/**
 * hs.c
 * CPSC 416 Assignment 3
 *
 * @author Mira Leung 
 * 
 * March 15, 2014
 *
 * Usage:
 * mpiexec -n <N PROCESSES> ./hs [ 1 ] for randomly assigned uids
 *
 * An implementation of Hirschberg-Sinclair's algorithm
 * for asynchronous ring leader election. Improvements on the HS algorithm are 
 * as follows (and marked in the code):
 * 1) Comparing to the maximum uid seen so far, instead of the process's own uid
 * 2) Reduce the number of election messages sent if a reply has already been received.
 *
 * Mesage complexity is in O(n log n):
 * * Worst case: less than 6n + 8n(ceiling{log n} - 1)
 * * * Proof: The minimum distance between two winners of a phase k-1 is 2^(k-1)+1, so the 
 *     total number of messags sent at phase k that is not the last phase is
 *     4(2^k * floor{n/(2^(k-1)+1)}) = 8n * floor{2^(k-1)/(2^(k-1)+1)} < 8n.
 *
 *     The total number of phases until the leader is elected is ceiling{log n} + 1, 
 *     including phase 0. 2n messages are sent in the last phase, and there are no reply messages.
 *     Hence, the total number of messages in the worst case is 
 *     4n + \sum^{\ceiling{log n} - 1}_{k=1}(4 * 2^k * n/(2^(k-1)+1)) + 2n
 *     <= 6n + 8n(\ceiling{log n} - 1). QED.
 * 
 * Time complexity: O(n). 
 * * The max total time for each phase k that is not the final one is 2(2^k), so the maximum
 *   total time required by phases 0 to the second-last one is 2(2^{ceiling(log n} + 1)), and
 *   the time for the final phase is n.
 *
 * Note: The message totals include the round of messages sent at the end to calculate the total.
 * The actual number sent during the algorithm is the final number minus n, where n = NUM.
 *
 * Total number of messages is roughly <= 6n + 8n * (ceiling{log n} - 1) + 3n (for 
 * the final round of message totals), which is 8n*(ceiling{log n}) + n \in O(n log n)
 *
 * The program checks if pnum is relatively coprime to and larger than size.
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include <fgmpi.h>


// Tags
#define TAG_ELECTION 2
#define TAG_REPLY 3
#define TAG_MSGNUM 4
#define TAG_DUMMY 5
#define TAG_IGNORE 6

#define SIZE_MSG 3

int ceiling_log2(unsigned long long x);
int gcd(int size, int pnum);

/** FG-MPI Boilerplate begins **/
int hs(int argc, char* argv[]);

FG_ProcessPtr_t binding_func(int argc, char** argv, int rank) {
  return (&hs);
}

FG_MapPtr_t map_lookup(int argc, char** argv, char* str) {
  return (&binding_func);
}

int main(int argc, char *argv[]) {
  FGmpiexec(&argc, &argv, &map_lookup);
  return 0;
}

/** FG-MPI Boilerplate ends **/


int hs(int argc, char *argv[]) {

  int lnum_sent = 0, lnum_recv = 0;
  int tnum_sent = 0, tnum_recv = 0;

  int rank, size;
  MPI_Init (&argc, &argv);  
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);  // get my pid
  MPI_Comm_size (MPI_COMM_WORLD, &size);  // get number of processes 
  MPI_Status status;
  MPI_Request request;
  
  int verbose = 0;
  if (argc == 2 && !strcmp(argv[1], "-v")) verbose = 1;

  int pnum = size * 1000000 + 1;

  int election_sendbuf[SIZE_MSG];
  srand(time(NULL) + rank);
  //int uid = ((rank+1)*pnum) % size;
  int uid = (rand() % pnum);
 
  int max_so_far = uid;
  int k = 0, d = 0;
  election_sendbuf[0] = uid, election_sendbuf[1] = k, election_sendbuf[2] = d;
  int left = rank-1;
  if (!rank) left = size-1;
  int right = (rank+1)%size;
  int recvbuf[SIZE_MSG];
  int left_recv_tag, right_recv_tag;
  int endLoopFlag = 0;
  int recvReplies[2][2] = {{0, 0}, {0, 0}}; // left, right; j, k

  int left_sendbuf[SIZE_MSG] = { max_so_far, k, d }, left_send_tag = TAG_ELECTION, left_send_dest = right;
  int right_sendbuf[SIZE_MSG] = { max_so_far, k, d }, right_send_tag = TAG_ELECTION, right_send_dest = left;

  int last = ceiling_log2((unsigned long long) size);

  MPI_Isend(election_sendbuf, SIZE_MSG, MPI_INT, left, TAG_ELECTION, MPI_COMM_WORLD, &request);
  MPI_Isend(election_sendbuf, SIZE_MSG, MPI_INT, right, TAG_ELECTION, MPI_COMM_WORLD, &request);
  lnum_sent+= 2;

  // Current leader is max_so_far
  while (k < last+1) {

    MPI_Recv(recvbuf, SIZE_MSG, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    lnum_recv++;
    k = recvbuf[1], d = recvbuf[2];
    if (status.MPI_TAG == TAG_IGNORE) {
      if (recvbuf[0] > max_so_far) max_so_far = recvbuf[0];
      break;
    }

    if (k > last) break;
 
    // Deal with messages received from the left
      if (status.MPI_SOURCE == left) {
    left_recv_tag = status.MPI_TAG;

      switch (left_recv_tag) {
        case TAG_ELECTION:    
          if (recvbuf[0] > uid) {
            if (recvbuf[2] <  (1 << k)) {
              left_sendbuf[2] = d + 1, left_sendbuf[0] = recvbuf[0], left_sendbuf[1] = k;
              left_send_tag = TAG_ELECTION, left_send_dest = right;
             } else if (d >= (1 << k)) {
              left_sendbuf[0] = recvbuf[0], left_sendbuf[1] = recvbuf[1];
              left_send_tag = TAG_REPLY, left_send_dest = left;
            } 
          } else if (recvbuf[0] == uid) {
            if (recvbuf[0] > max_so_far) max_so_far = recvbuf[0];
            left_sendbuf[0] = uid, left_sendbuf[1] = k+1, left_sendbuf[2] = 1;
            left_send_tag = TAG_ELECTION, left_send_dest = right;
          } else {
            break; 
          }
          lnum_sent++;
          MPI_Isend(left_sendbuf, SIZE_MSG, MPI_INT, left_send_dest, left_send_tag, MPI_COMM_WORLD, &request);
          break;

      case TAG_REPLY:
          if (recvbuf[0] != max_so_far) { // Improvement #1: compare to max_so_far instead of own uid
            if (recvbuf[0] > max_so_far) max_so_far = recvbuf[0];
            left_sendbuf[0] = recvbuf[0], left_sendbuf[1] = recvbuf[1];
            left_send_dest = right, left_send_tag = TAG_REPLY;
            recvReplies[0][0] = recvbuf[0], recvReplies[0][1] = recvbuf[1];
          } else {
            if (recvbuf[0] > max_so_far) max_so_far = recvbuf[0];
            left_sendbuf[0] = recvbuf[0], left_sendbuf[1] = k + 1, left_sendbuf[2] = 1;
            lnum_sent+=2;
           MPI_Isend(left_sendbuf, SIZE_MSG, MPI_INT, left, TAG_ELECTION, MPI_COMM_WORLD, &request); 
           MPI_Isend(left_sendbuf, SIZE_MSG, MPI_INT, right, TAG_ELECTION, MPI_COMM_WORLD, &request);
          }
          break;

      case TAG_IGNORE: break;
      default: endLoopFlag = 1; break;
    }
      }
   

    // Deal with right received values
      else if (status.MPI_SOURCE == right) {
        right_recv_tag = status.MPI_TAG;
    switch (right_recv_tag) {
      case TAG_ELECTION:
          if (recvbuf[0] > uid) { 
            if (recvbuf[0] > max_so_far) max_so_far = recvbuf[0];
            if (d < (1 << k)) {
              right_sendbuf[2] = recvbuf[2]+1, right_sendbuf[0] = recvbuf[0], right_sendbuf[1] = recvbuf[1];
              right_send_tag = TAG_ELECTION, right_send_dest = left;
            } else if (d >= (1 << k)) {
              right_sendbuf[0] = recvbuf[0], right_sendbuf[1] = recvbuf[1];
              right_send_tag = TAG_REPLY, right_send_dest = right;
            }
          } else if (recvbuf[0] == uid) {
            if (recvbuf[0] > max_so_far) max_so_far = recvbuf[0];
            right_sendbuf[0] = uid, right_sendbuf[1] = k+1, right_sendbuf[2] = 1;
            right_send_tag = TAG_ELECTION, right_send_dest = left;
            if (k >= last && recvbuf[0] == max_so_far) endLoopFlag = 1;
          } else {
            break;           
          }
          lnum_sent++;
          MPI_Isend(right_sendbuf, SIZE_MSG, MPI_INT, right_send_dest, right_send_tag, MPI_COMM_WORLD, &request);
          break;

    case TAG_REPLY:
        if (recvbuf[0] != max_so_far) {
          if (recvbuf[0] > max_so_far) max_so_far = recvbuf[0];
          right_sendbuf[0] = recvbuf[0], right_sendbuf[1] = recvbuf[1];
          right_send_tag = TAG_REPLY, right_send_dest = left;
           recvReplies[1][0] = recvbuf[0], recvReplies[1][1] = recvbuf[1];

        } else {
          if (recvbuf[0] > max_so_far) max_so_far = recvbuf[0];
            if (recvReplies[1][0] == recvbuf[0] && recvReplies[1][1] == recvbuf[1]) {
          left_sendbuf[0] = uid,  
           left_sendbuf[1] = k+1, right_sendbuf[2] = left_sendbuf[2] = 1;
          lnum_sent+=2;
          MPI_Isend(left_sendbuf, SIZE_MSG, MPI_INT, left, TAG_ELECTION, MPI_COMM_WORLD, &request);
          MPI_Isend(left_sendbuf, SIZE_MSG, MPI_INT, right, TAG_ELECTION, MPI_COMM_WORLD, &request);
           } else {
            recvReplies[1][0] = recvbuf[0], recvReplies[1][1] = recvbuf[1];
          }
        break;
      }

    case TAG_IGNORE: break;
    default: endLoopFlag = 1; break;
    }
  }
      
    if (endLoopFlag) {
      break;
    }
 }
  
  int msgBuf[3] = {max_so_far, 0, 0}, msgRecv[3];

  // Election is over - tell the other processes
  MPI_Isend(msgBuf, SIZE_MSG, MPI_INT, left, TAG_IGNORE, MPI_COMM_WORLD, &request);
  MPI_Isend(msgBuf, SIZE_MSG, MPI_INT, right, TAG_IGNORE, MPI_COMM_WORLD, &request);
  MPI_Request_free(&request);

  if (verbose) printf("rank=%d, id=%d, leader=%d, mrcvd=%d, msent=%d\n", rank, uid, (max_so_far == uid), lnum_recv, lnum_sent);
  if (max_so_far == uid) {
      msgBuf[0] = lnum_recv, msgBuf[1] = lnum_sent, msgBuf[2] = uid;
      max_so_far = uid;
      MPI_Isend(msgBuf, SIZE_MSG, MPI_INT, right, TAG_MSGNUM, MPI_COMM_WORLD, &request); 
  } else {
     // Non-leaders, do a receive and a send for their message totals as well
      MPI_Recv(msgRecv, SIZE_MSG, MPI_INT, left, TAG_MSGNUM, MPI_COMM_WORLD, &status);
      // increase its count by 1 for each receive and send
      msgBuf[0] = msgRecv[0]+ lnum_recv, msgBuf[1] = msgRecv[1]+ lnum_sent;
      max_so_far = msgBuf[2];
      MPI_Isend(msgBuf, SIZE_MSG, MPI_INT, right, TAG_MSGNUM, MPI_COMM_WORLD, &request); 
  }

  // Leader receives/prints total number of messages received and sent
  if (max_so_far == uid) {
    MPI_Recv(msgRecv, 3, MPI_INT, left, TAG_MSGNUM, MPI_COMM_WORLD, &status); 
    tnum_recv = msgRecv[0];
    tnum_sent = msgRecv[1];
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

