#ifndef HMR_UI_H
#define HMR_UI_H

/* Print the text with the current time at the beginning */
void time_print(const char *str, ...);

/* Exit the program with an error code, and print an error info */
void time_error(int exitCode, const char *str, ...);

#endif // HMR_UI_H
