./floattest <./testes/caso-00.in >saida0.out
./floattest <./testes/caso-01.in >saida1.out
./floattest <./testes/caso-02.in >saida2.out
./floattest <./testes/caso-03.in >saida3.out
./floattest <./testes/caso-04.in >saida4.out
./floattest <./testes/caso-05.in >saida5.out
./floattest <./testes/caso-06.in >saida6.out
./floattest <./testes/caso-07.in >saida7.out
./floattest <./testes/caso-08.in >saida8.out
./floattest <./testes/caso-09.in >saida9.out
echo "" > resp.txt
echo "------------------------[ CASO 0 ]-----------------------------" >> resp.txt
echo ""  >> resp.txt
diff -w ./testes/caso-00.out saida0.out >> resp.txt
echo "" >> resp.txt
echo "------------------------[ CASO 1 ]-----------------------------" >> resp.txt
echo "" >> resp.txt
diff -w ./testes/caso-01.out saida1.out  >> resp.txt
echo "" >> resp.txt
echo "------------------------[ CASO 2 ]-----------------------------" >> resp.txt
echo "" >> resp.txt
diff -w ./testes/caso-02.out saida2.out  >> resp.txt
echo "" >> resp.txt
echo "------------------------[ CASO 3 ]-----------------------------" >> resp.txt
echo "" >> resp.txt
diff -w ./testes/caso-03.out saida3.out  >> resp.txt
echo "" >> resp.txt
echo "------------------------[ CASO 4 ]-----------------------------" >> resp.txt
echo "" >> resp.txt
diff -w ./testes/caso-04.out saida4.out  >> resp.txt
echo "" >> resp.txt
echo "------------------------[ CASO 5 ]-----------------------------" >> resp.txt
echo "" >> resp.txt
diff -w ./testes/caso-05.out saida5.out  >> resp.txt
echo "" >> resp.txt
echo "------------------------[ CASO 6 ]-----------------------------" >> resp.txt
echo "" >> resp.txt
diff -w ./testes/caso-06.out saida6.out  >> resp.txt
echo "" >> resp.txt
echo "------------------------[ CASO 7 ]-----------------------------" >> resp.txt
echo "" >> resp.txt
diff -w ./testes/caso-07.out saida7.out  >> resp.txt
echo "" >> resp.txt
echo "------------------------[ CASO 8 ]-----------------------------" >> resp.txt
echo "" >> resp.txt
diff -w ./testes/caso-08.out saida8.out  >> resp.txt
echo "" >> resp.txt
echo "------------------------[ CASO 9 ]-----------------------------" >> resp.txt
echo "" >> resp.txt
diff -w ./testes/caso-09.out saida9.out >> resp.txt
echo "" >> resp.txt
cat resp.txt
