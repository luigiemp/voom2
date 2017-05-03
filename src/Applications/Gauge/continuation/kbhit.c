/**
 * Unix C implementation of kbhit(), based on
 * Morgan McGuire, morgan@cs.brown.edu
 * <http://www.flipcode.org/cgi-bin/fcarticles.cgi?show=64166>
 *
 * gcc -c kbhit.c -o kbhit.o
 */

/* --- self-identity --- */
#include "kbhit.h"

/* fileno setbuf stdin */
#include <stdio.h>

/* NULL */
#include <stddef.h>

/* termios tcsetattr tcgetattr TCSANOW */
#include <termios.h>

/* ioctl FIONREAD ICANON ECHO */
#include <sys/ioctl.h>

static int initialized = 0;
static struct termios original_tty;


int kbhit() 
{
  if(!initialized)
  {
    kbinit();
  }

  int bytesWaiting;
  ioctl(fileno(stdin), FIONREAD, &bytesWaiting);
  return bytesWaiting;
}

/* Call this just when main() does its initialization. */
/* Note: kbhit will call this if it hasn't been done yet. */
void kbinit()
{
  struct termios tty;
  tcgetattr(fileno(stdin), &original_tty);
  tty = original_tty;

  /* Disable ICANON line buffering, and ECHO. */
  tty.c_lflag &= ~ICANON;
  //tty.c_lflag &= ~ECHO;
  tcsetattr(fileno(stdin), TCSANOW, &tty);

  /* Decouple the FILE*'s internal buffer. */
  /* Rely on the OS buffer, probably 8192 bytes. */
  setbuf(stdin, NULL);
  initialized = 1;
}

/* Call this just before main() quits, to restore TTY settings! */
void kbfini()
{
  if(initialized)
  {
    tcsetattr(fileno(stdin), TCSANOW, &original_tty);
    initialized = 0;
  }
}


