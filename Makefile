default:
	@echo "Usage:"
	@echo "    all: compiles the program"
	@echo "    lib: compiles pfft library"
	@echo "    clean: removes compilated objects & executable files"

all:
	g++ pfft.h pfft.cpp D_test.cpp -o dt -Wall -Wextra -g
	g++ pfft.h pfft.cpp D2_test.cpp -o it -Wall -Wextra -g
	g++ pfft.h pfft.cpp fft_test.cpp -o pftest -Wall -Wextra -g

lib:
	g++ pfft.h pfft.cpp -Wall -Wextra

clean:
	rm pftest dt
