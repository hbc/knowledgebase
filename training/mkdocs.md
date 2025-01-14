# MK-Docs

Basic prerequisites
	-Python
	-Github CLI
	-VS Code


Download Visual Studio Code

https://code.visualstudio.com/download

1. Open new terminal within VS-Code
2. Navigate to where you want to host the GitHub repo.
3. Create a new GitHub repo you want to work on or clone an existing one: git clone link_to_repo.git
4. Cd to the GitHub repo folder
5. Lets make a virtual python environment: python -m venv venv 
6. Source the python environment you just created: source venv/bin/activate
7. Need pip installed check: pip —version
8. Install mkdocs: pip install mkdocs-material
9. Open visual code from here: code .
10. Open terminal within vscode
11. To open up the website mkdocs serve
12.  To change to the “material theme” open mkdocs.yml file and below site_name: My Docs, type 
		site_name: My Docs
		theme:
			name: material
13. Save it
14.  To deploy: type mkdocs serve on terminal. It will restart in the same host.
15. We can change the appearance  of the site by editing and adding plugins to mkdocs.yml (google for common settings)
16. To add another page. Go inside docs add another nage: eg page2.md
17. Add 
	# Page 2
	
	## sub heading

	Text inside
18. Save it
19. In the terminal type: git add .
20. Git commit -m $’updating instructor’
	git push origin main
	git config http.postBuffer 524288000 #(if error comes up related to html)
	git pull && git push
