#ifndef API_H
#define API_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <curl/curl.h>

typedef struct {
  char *memory;
  size_t size;
} memory;

size_t write_memory(void *contents, size_t size, size_t nmemb, void *userp);
void get_vectors(CURL *curl, memory *chunk, double *vectors, int year, const char *id);

#endif
