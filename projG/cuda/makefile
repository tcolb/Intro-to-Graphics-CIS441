CC = nvcc
CFLAGS = `pkg-config --cflags --libs opencv`

project1G: project1G.cu
	$(CC) -o $@ $< $(CFLAGS)

clean:
	rm project1G
