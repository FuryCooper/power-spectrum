#ifndef IO_H_
#define IO_H_

void read_param();
#ifdef NEUTRINO
void load_neutrino_ratio();
#endif // NEUTRINO
void load_snapshot(int rep);
void output();

#endif // !IO_H_
