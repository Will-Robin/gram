ROOT=../gram

all: index.html

index.html: $(GRAM)
	pdoc ../gram -d "numpy" -o ./

clean:
	rm -rf gram
	rm gram.html
	rm index.html
	rm search.js
