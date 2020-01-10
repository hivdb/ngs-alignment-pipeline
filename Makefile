build:
	@docker build . -t hivdb/five-prime-alignment-exp:latest

shell:
	@docker run -it --rm -v $(PWD):/workspace hivdb/five-prime-alignment-exp:latest bash

release:
	@docker push hivdb/five-prime-alignment-exp:latest
