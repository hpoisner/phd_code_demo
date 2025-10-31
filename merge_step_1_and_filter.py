import subprocess
import argparse

parser = argparse.ArgumentParser(description = 'Arguments for generating commands file')

parser.add_argument('--output_folder', '-of', type=str, required=True, help='Path to output folder, must use double quotes')
parser.add_argument('--output_file', '-ofi', type=str, required=True, help='Command File name, must use double quotes')

args = parser.parse_args()

OUTPUT_FOLDER = args.output_folder
outfile = open(args.output_file, "w")

    
# Construct dx run command
cmd = f'''dx run app-swiss-army-knife \
--priority normal \
--instance-type "mem1_ssd1_v2_x36" \
--tag="merge_gwas_step1" \
-iin=file-J0J3zKjJvjB83k4QBB1f72Zp \
-iin=file-J0J42XjJVy21xQ5G85ZJy20k \
-iin=file-J0J413jJgqY4fkKV0g3YYJk6 \
-iin=file-J0J3yG0JPqY83k4QBB1f72FB \
-iin=file-J0J3qv0JZJbv23P6Xp7yY6pj \
-iin=file-J0J3y2QJJgg1ZFYqKFK7Yy42 \
-iin=file-J0J3yFQJygVqQG2qX9X82XKg \
-iin=file-J0J3xq8JzQkQ7J45bF3X5Zvf \
-iin=file-J0J41XjJ5qv2kp9ZK5Qxb7xv \
-iin=file-J0J3v68JqjX2zKK5qQjv8Z81 \
-iin=file-J0J43XQJqgVjfB1G9zBy1PZz \
-iin=file-J0J3yX8JFQbf0Pj5fXX4P6v0 \
-iin=file-J0J3vXQJvBPG23P6Xp7yY6x9 \
-iin=file-J0J3P70J2g2zj7F7f1x7PP24 \
-iin=file-J0J41v0J52J7KjPPbq1j7QPp \
-iin=file-J0J408QJ4Bq61ZZ5g05X2VYk \
-iin=file-J0J404QJ4qj37Bb1X3JVKpxP \
-iin=file-J0J3zG0Jyfvf4kQQvXxB14v5 \
-iin=file-J0J3z80JVQ5VfPJB6YvGfvvq \
-iin=file-J0J41kjJBKxqKg978BP1Zkvy \
-iin=file-J0J3zbjJ4q4f1ZZ5g05X2VQQ \
-iin=file-J0J41X8JygX83k4QBB1f72z7 \
-iin=file-J0J44F0JG1190JxG9QFq41k1 \
-iin=file-J0J3zKjJvjBB6ZbY28KgP2vP \
-iin=file-J0J42XjJVy2G23P6Xp7yY7zk \
-iin=file-J0J413jJgqYJk8PBKxypBg3f \
-iin=file-J0J3yG0JPqY64kQQvXxB14gZ \
-iin=file-J0J3qv0JZJbxx2z8x972XyZ9 \
-iin=file-J0J3y2QJJggKz020kKYxk72b \
-iin=file-J0J3yFQJygVj3k4QBB1f72F0 \
-iin=file-J0J3xq8JzQkgxYg328vX85zY \
-iin=file-J0J41XjJ5qv0FKf2j7V1VgyF \
-iin=file-J0J3v68JqjXBzQfyVq0PzB9J \
-iin=file-J0J43XQJqgVbkFxG926g3Z54 \
-iin=file-J0J3yX8JFQbf1ZZ5g05X2VBV \
-iin=file-J0J3vXQJvBPKpk98J8BPkg4P \
-iin=file-J0J3P70J2g2z7fpFXY012vPk \
-iin=file-J0J41v0J52J0G3B3345xFX2K \
-iin=file-J0J408QJ4Bq6pF3k8yG7K1zX \
-iin=file-J0J404QJ4qjP1QVK0B0V3Vv1 \
-iin=file-J0J3zG0Jyfvzq0b9zp3Yf7XG \
-iin=file-J0J3z80JVQ5yQ0Fq8P2FP9VP \
-iin=file-J0J41kjJBKxzX79Z9v5JQ46P \
-iin=file-J0J3zbjJ4q4QQ9JGygqVVBGG \
-iin=file-J0J41X8JygX7KjPPbq1j7QK4 \
-iin=file-J0J44F0JG116yqp9X0xBkZ23 \
-iin=file-J0J3zKjJvjBPZG0Fvj9gQyK4 \
-iin=file-J0J42XjJVy28fB1G9zBy1BZ4 \
-iin=file-J0J413jJgqY8VP8XpZ7fJFK6 \
-iin=file-J0J3yG0JPqY3Xvk4865FXGJb \
-iin=file-J0J3qv0JZJbY7Bb1X3JVGzG4 \
-iin=file-J0J3y2QJJgg64kQQvXxB14Z4 \
-iin=file-J0J3yFQJygVjVgQqkv81p0P8 \
-iin=file-J0J3xq8JzQkxxG1KgG3k7Zgx \
-iin=file-J0J41XjJ5qv8fB1G9zBy1BGZ \
-iin=file-J0J3v68JqjX9bJXXgKZq3jVK \
-iin=file-J0J43XQJqgVyzYy6k6f9jKBQ \
-iin=file-J0J3yX8JFQbf0kq029x5yZyP \
-iin=file-J0J3vXQJvBP1Zgq6pzVKBj5j \
-iin=file-J0J3P70J2g2Vxyk0K5JY8y79 \
-iin=file-J0J41v8J52JGZ8jyVV20F1y5 \
-iin=file-J0J408QJ4BqP1QVK0B0V3VvZ \
-iin=file-J0J404QJ4qj1fPJB6YvGfx1b \
-iin=file-J0J3zG0Jyfvv2kkGVzP8V3q6 \
-iin=file-J0J3z80JVQ5ZfkKV0g3YYJPZ \
-iin=file-J0J41kjJBKxfyqp9X0xBkYg9 \
-iin=file-J0J3zbjJ4q4zj7F7f1x7Qj4J \
-iin=file-J0J41X8JygXG6yv9vxB3yJzV \
-iin=file-J0J44F0JG11Px2Y0PvvxgkzB \
-iin=file-J0JjVq8J2p5Y2P94v3fQ2Q5k \
-icmd="plink --merge-list mergelist.txt --make-bed --thin-count 250000 --out ukbb_wgs_merged" \
--destination "{OUTPUT_FOLDER}" \
-y --brief \
--ignore-reuse\n'''

# Write command to output file
outfile.write(cmd)
