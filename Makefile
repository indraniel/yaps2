.PHONY: clean

init:
	pip install -r requirements.txt

test:
	echo "Please implement me!"

clean:
	find ./yaps2 -name "*.pyc" -exec rm {} \;
