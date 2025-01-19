---
layout: page
title: Programa
subtitle: Introdução à Programação de Computadores para Biologia
---


## 2o semestre de 2016

**Horário:**    Sexta-Feira – 14:00 às 16:00 horas  
**Professora:** Tatiana Teixeira Torres, e-mail: tttorres(bio)ib.usp.br

OBS: substituir (bio) por @ 

## Programa:

![Programa](/img/programa2018.png)

## Formas de avaliação:

1- Exercícios semanais (E): a cada novo tema de aula, será disponibilizado um exercício que deverá ser entregue no prazo de sete dias. A entrega dos exercícios será realizada via e-mail (endereço acima). O atraso na entrega implica na redução da nota pela metade. Serão 10 exercícios ao longo do semestre. 

2- Projeto (P): a cada grupo de três alunos, será atribuído um objetivo de pesquisa. O projeto deverá descrever como os alunos pretender abordar o problema utilizando conhecimento de algoritmos e programação de computadores. O projeto deverá ser estruturado em Introdução e justificativa, Objetivos e Métodos e não deverá ultrapassar quatro páginas.

3- Manuscrito (M): o manuscrito descreverá os resultados obtidos no desenvolvimento do projeto. Juntamente com o manuscrito, deverão ser entregues os scripts preparados pelos alunos. Cada aluno deverá preparar pelo menos um script. Os alunos serão avaliados em grupo pelo manuscrito (peso 4) e individualmente pelos scripts (peso seis).

## Critério para aprovação e pesos das avaliações:

{% highlight perl linenos %}  
#!/usr/bin/perl

my ($E, $P, $M, $nota_final);
($E, $P, $M) = @ARGV;

$nota_final = ((2*$E)+(3*$P)+(5*$M))/10;

if ($nota_final >= 5) {
  print "Aluno aprovado\n";
} else {
  print "Aluno reprovado\n";
}

#conversao de nota para conceitos, apenas para a Pos-graduacao
if ($nota_final < 5) {
  print "R, Reprovado, sem direito a credito\n";
} elsif ($nota_final <= 7.0) {
  print "C, Regular, com direito a credito\n";
} elsif ($nota_final <= 8.5) {
  print "B, Bom, com direito a credito\n";
} else {
  print "A, Excelente, com direito a credito\n";
}
exit;

{% endhighlight %}


## Bibliografia:

- Curtis Jamison, D. “Perl Programming for Biologists”, John Wiley & Sons, Inc., 2003.
- Wall, L.; Christiansen, T.; Orwant, J. “Programming Perl”, O'Reilly Media, 2000.
- Tisdall, J. “Beginning Perl for Bioinformatics”, O'Reilly Media, 2001.
- Setubal J. C.; Meidanis J. "Introduction to Computational Molecular Biology", Brooks/Cole Pub Co, 1997.
