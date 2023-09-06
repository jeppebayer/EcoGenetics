<!DOCTYPE html>

<html>

<!--    NAVIGATION BAR-->
<?php include "header.html";?>
<?php include "title.html";?>

<?php
    $aResult = array();
    if( !isset($_POST['code']) ) { $aResult['error'] = 'ERROR: No code passed to input_validation.php';}
    $code=$_POST["code"];
    if( !isset($_POST['kmer_length']) ) { echo shell_exec('echo ERROR: No kmer_length passed to input_validation.php >> user_data/$code/input_validation.log');}
    if( !isset($_POST['read_length']) ) { echo shell_exec('echo ERROR: No read_length passed to input_validation.php >> user_data/$code/input_validation.log');}
    if( !isset($_POST['description']) ) { echo shell_exec('echo ERROR: No description passed to input_validation.php >> user_data/$code/input_validation.log');}

    if( !isset($_POST['max_kmer_cov']) ) { echo shell_exec('echo ERROR: No max_kmer_cov passed to input_validation.php >> user_data/$code/input_validation.log');}

    $description = $_POST["description"];
    $kmer_length = $_POST["kmer_length"];
    $read_length = $_POST["read_length"];
    $max_kmer_cov = $_POST["max_kmer_cov"];


    $url="analysis.php?code=$code";
    $run_url="run_now.php";
    $filename="user_uploads/$code";
    
    $back_button= "<form action=\"./\" method=GET><button type=\"submit\" class=\"center btn btn-danger\">Back</button></form>";
    //$continue_button= "<form action=\"$url\"><input type=\"hidden\" name = \"code\" value=\"$code\"><button type=\"submit\" class=\"center btn btn-success\">Continue</button></form>";
    
    $continue_button= "<form action=\"$run_url\" method=\"post\"><input type=\"hidden\" name = \"code\" value=\"$code\">   <input type=\"hidden\" name=\"kmer_length\" value=\"$kmer_length\"> <input type=\"hidden\" name=\"description\" value=\"$description\"> <input type=\"hidden\" name=\"max_kmer_cov\" value=\"$max_kmer_cov\">  <input type=\"hidden\" name=\"read_length\" value=\"$read_length\">   <button type=\"submit\" class=\"center btn btn-success\">Continue</button></form>";
    
    
    if (!file_exists ($filename)) {
        echo "<div class=\"alert center alert-danger\" role=\"alert\">No file uploaded</div>";
        echo "$back_button";
        exit;
    }
    
    $myfile = fopen($filename, "r") or die("Unable to open file!");
    
    $line_counter=0;
    $num_columns=0;
    $bin_history=array();
    $consistent=true;
    while(!feof($myfile)) {
        $bin_counter=0;
        $line=fgets($myfile);
        // $line =trim(preg_replace( '/\s+/', ' ', $line ));
        if ($line=="" or $line==" ") {
            continue;
        }
        $array=array_map("trim",explode(' ',$line));
        //var_dump($array);
        $bin_counter=count($array);
        //echo $bin_counter . "<br>";
        //echo $num_columns. "<br>";
        
        if ($num_columns != 0 and $num_columns != $bin_counter) {
            $consistent=false;
            // echo $line;
            // var_dump($array);
        }
        $num_columns=$bin_counter;
        $line_counter=$line_counter+1;
        $bin_history[]=$bin_counter;
    }
    fclose($myfile);
    
    if ($consistent) {
        if ($line_counter >= 50 and $num_columns == 2) {
            echo "<div class=\"alert center alert-success\" role=\"alert\">Great! File was uploaded and has acceptable dimensions: $line_counter rows by $num_columns columns</div>";
            echo "<div style=\"margin-left:1%;\"><div class=\"col-sm-1\">";
            echo "$back_button";
            echo "</div><div class=\"col-sm-1\">";
            echo "$continue_button";
            echo "</div></div>";
        } else {
            if ($num_columns != 2) {
                echo "<div class=\"alert center alert-danger\" role=\"alert\">File was uploaded but it has $num_columns column(s). The file must have 2 columns separated by a single space, which is the default in Jellyfish</div>";    
            }
            if ($line_counter < 50) {
                echo "<div class=\"alert center alert-danger\" role=\"alert\">File was uploaded but it only has $line_counter rows, are you sure this is the right file?</div>";        
            }
            echo "$back_button";
        }
    }
    else {
        echo "<div class=\"alert center alert-danger\" role=\"alert\">All lines in file must have the same number of elements (separated by single spaces). This file had the following numbers of elements per line: ";
        foreach ($bin_history as $num)
            echo $num . ", ";
        echo "</div>";
        echo "$back_button";
    }
   
    
    
    
?>
</body>
</html>
