CC		= gcc
CFLAGS	= -Wall -Wextra -Werror -std=c99
SRC		= $(wildcard src/*.c)
OBJ     = $(SRC:.c=.o)
TARGET	= n-body

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -lm -lcurl -o $@ $^

src/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@

run: $(TARGET)
	./$(TARGET)

clean:
	rm $(TARGET) src/*.o *.json

.PHONY: clean run'
